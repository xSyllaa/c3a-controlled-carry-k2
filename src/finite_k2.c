/* finite_k2.c — k=2 min-cost residue histogram for the controlled-carry construction.
   Mask M = <1971,2016,2100,2628,2688,2800> ∩ [0,26972], base q=27022, residues mod q^2.
   c2(r) = min over lifts e = d0 + d1*q ≡ r (mod q^2), d_i ∈ M-M, of kappa(d0)+kappa(d1),
   kappa(d) = min{a+b : a,b ∈ M, a-b = d}.  The histogram does NOT depend on x.
   Usage: ./finite_k2 OUT_HIST.tsv OUT_META.json   (needs ~1.5 GB RAM)
   The meta file records x-independent facts plus binary64 Z2/P+ at the CERTIFICATE x (info only;
   certified enclosures are in certificates/cert_primary_ell_1e18.json).
   Credit: Sylla / xSyllaa. (2026-10-19-c3a-xopt: fixed double-fclose, output paths as args;
   hostile 119281555 polish: info x switched from retired 0x1.ffde827adc0fep-1 to x*, meta made x-explicit.) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

enum { B=26972, QBASE=27022, INF=1000000000 };
static const int GENS[] = {1971,2016,2100,2628,2688,2800};
static const int NGENS = 6;
static int *cost; static unsigned char *inmask;
static int *Dvals, *Dcosts, nD;
static double Pplus_cache;
#define XINFO_HEX "0x1.ffde1bafccc06p-1"   /* certificate x* (info only; histogram is x-independent) */

static void build(void) {
  inmask = calloc(B+1,1); inmask[0]=1;
  for(int n=0;n<=B;n++) if(inmask[n]) for(int i=0;i<NGENS;i++){
    int g=GENS[i]; if(n+g<=B) inmask[n+g]=1;
  }
  cost = malloc((B+1)*sizeof(int));
  for(int d=0;d<=B;d++) cost[d]=INF; cost[0]=0;
  int *M=malloc((B+1)*sizeof(int)); int nM=0;
  for(int n=0;n<=B;n++) if(inmask[n]) M[nM++]=n;
  for(int i=0;i<nM;i++) for(int j=0;j<nM;j++){
    int a=M[i],b=M[j],d=a-b; if(d<0)d=-d;
    int c=a+b; if(c<cost[d]) cost[d]=c;
  }
  Dvals=malloc((2*(B+1)+2)*sizeof(int));
  Dcosts=malloc((2*(B+1)+2)*sizeof(int));
  nD=0; Dvals[nD]=0; Dcosts[nD]=0; nD++;
  for(int d=1;d<=B;d++) if(cost[d]<INF){
    Dvals[nD]=d; Dcosts[nD]=cost[d]; nD++;
    Dvals[nD]=-d; Dcosts[nD]=cost[d]; nD++;
  }
  unsigned char *sum=calloc(2*B+1,1);
  for(int i=0;i<nM;i++) for(int j=0;j<nM;j++) sum[M[i]+M[j]]=1;
  free(M);
  double x = strtod(XINFO_HEX, 0);
  double Pp=0; for(int s=0;s<=2*B;s++) if(sum[s]) Pp+=pow(x,s);
  Pplus_cache=Pp; free(sum);
  fprintf(stderr,"|D|=%d P+=%.17g\n", nD, Pp);
}

/* Build c2 array and histogram m[c] = number of residues with cost c.
   Also track cost of 0 and of h=Q2/2. */
static uint16_t *build_c2(int *hist, int maxc_cap, int *c0_out, int *ch_out) {
  long Q2 = (long)QBASE*(long)QBASE;
  uint16_t *c2 = malloc(Q2 * sizeof(uint16_t));
  if(!c2){perror("malloc"); exit(1);}
  memset(c2, 0xFF, Q2*sizeof(uint16_t));
  for(int i=0;i<nD;i++){
    int d0=Dvals[i], c0=Dcosts[i];
    long base = d0; if(base<0) base += Q2;
    for(int j=0;j<nD;j++){
      int d1=Dvals[j], c1=Dcosts[j];
      long val = base + (long)d1 * QBASE;
      long r = val % Q2; if(r<0) r+=Q2;
      int cst = c0+c1;
      if(cst>65534) continue;
      if((uint16_t)cst < c2[r]) c2[r]=(uint16_t)cst;
    }
    if((i&4095)==0) fprintf(stderr,"\r build i=%d/%d", i, nD);
  }
  fprintf(stderr,"\nhistogram...\n");
  memset(hist, 0, maxc_cap*sizeof(int));
  for(long r=0;r<Q2;r++){
    if(c2[r]!=0xFFFF){
      int c=c2[r];
      if(c>=maxc_cap){fprintf(stderr,"cost overflow %d\n",c); exit(1);}
      hist[c]++;
    }
  }
  *c0_out = c2[0];
  long h = Q2/2;
  *ch_out = (c2[h]==0xFFFF)? -1 : (int)c2[h];
  return c2;
}

int main(int argc, char **argv){
  if(argc<3){fprintf(stderr,"usage: %s hist.tsv meta.json\n",argv[0]); return 2;}
  double x = strtod(XINFO_HEX, 0);
  build();
  int maxc_cap = 70000;
  int *hist = calloc(maxc_cap, sizeof(int));
  int c0, ch;
  uint16_t *c2 = build_c2(hist, maxc_cap, &c0, &ch);
  long Q2 = (long)QBASE*(long)QBASE;
  double Z=0; int ncosts=0;
  for(int c=0;c<maxc_cap;c++) if(hist[c]){ Z += hist[c]*pow(x,c); ncosts++; }
  fprintf(stderr,"Z2=%.17g distinct_costs=%d c0=%d ch=%d\n", Z, ncosts, c0, ch);

  /* Save histogram for Python finite-ell */
  FILE *f = fopen(argv[1],"w");
  fprintf(f, "# cost\tmultiplicity\n");
  for(int c=0;c<maxc_cap;c++) if(hist[c]) fprintf(f, "%d\t%d\n", c, hist[c]);
  fclose(f);
  int maxc=0; long covered=0; for(int c=0;c<maxc_cap;c++) if(hist[c]){ maxc=c; covered+=hist[c]; }
  FILE *m = fopen(argv[2],"w");
  fprintf(m, "{\"B\":%d,\"q\":%d,\"k\":2,\"Q2\":%ld,\"covered\":%ld,\"c0\":%d,\"ch\":%d,\"distinct_costs\":%d,\"max_cost\":%d,"
             "\"x_hex_certificate\":\"%s\",\"Z2_at_x_binary64_info\":%.17g,\"Pplus_at_x_binary64_info\":%.17g,"
             "\"note\":\"histogram is x-independent; Z2/Pplus here are binary64 info only, certified enclosures are in cert_primary_ell_1e18.json\"}\n",
    B, QBASE, Q2, covered, c0, ch, ncosts, maxc, XINFO_HEX, Z, Pplus_cache);
  fclose(m);
  printf("wrote hist + meta (x-independent histogram; Z2/P+ info at certificate x %s, binary64) Z2=%.17g sqrtZ2=%.17g theta_inf_info=%.17g\n", XINFO_HEX, Z, sqrt(Z), 1+log(sqrt(Z)/Pplus_cache)/log(QBASE));
  free(c2); free(hist);
  return 0;
}
