/* k3_auto_A.c — Builder A (independent of B): carry-automaton k=3 min-cost histogram mod q^3.
   r = r0 + r1 q + r2 q^2; lift digits d_i = r_i - g_i + q g_{i+1}, carries g_i in {0,-1}, g_0=0, g_3 free;
   cost = sum kappa(d_i). A_g(r0,r1) = min_{g1} kappa(r0+q g1)+kappa(r1-g1+q g); T_g(r2) = min_{g3} kappa(r2-g+q g3);
   c3 = min_g A_g + T_g.  Also emits k=2 histogram min_g A_g (must equal published k2_cost_hist.tsv).
   No c2 array, no D x D scatter (unlike builder B).  Usage: k3_auto_A OUT_K3.tsv OUT_K2.tsv */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#define BM 26972
#define QQ 27022
#define INF (1L<<30)
#define CM 200000
static long kap[2*BM+1];
static long K(long d){ return (d<-BM||d>BM)?INF:kap[d+BM]; }
static long e[2][2][QQ];
static int cmp(const void*a,const void*b){ uint32_t x=*(const uint32_t*)a,y=*(const uint32_t*)b; return (x>y)-(x<y); }
static long mn(long a,long b){return a<b?a:b;}
int main(int argc,char**argv){
  static char in[BM+1]; in[0]=1; int g[6]={1971,2016,2100,2628,2688,2800};
  for(int n=1;n<=BM;n++) for(int i=0;i<6;i++) if(n>=g[i]&&in[n-g[i]]){in[n]=1;break;}
  static int M[BM+1]; int nM=0; for(int n=0;n<=BM;n++) if(in[n]) M[nM++]=n;
  for(int i=0;i<=2*BM;i++) kap[i]=INF;
  for(int i=0;i<nM;i++) for(int j=0;j<nM;j++){ int d=M[i]-M[j]+BM; if(M[i]+M[j]<kap[d]) kap[d]=M[i]+M[j]; }
  for(int i=0;i<2;i++) for(int o=0;o<2;o++) for(int r=0;r<QQ;r++) e[i][o][r]=K((long)r+i-(long)QQ*o);
  long Q2=(long)QQ*QQ;
  uint32_t *key=malloc(Q2*4); uint64_t *h2=calloc(CM,8);
  #pragma omp parallel
  { uint64_t *loc=calloc(CM,8);
    #pragma omp for schedule(static)
    for(long r1=0;r1<QQ;r1++) for(long r0=0;r0<QQ;r0++){
      long A[2];
      for(int o=0;o<2;o++){ long best=INF; for(int c1=0;c1<2;c1++) best=mn(best,e[0][c1][r0]+e[c1][o][r1]); A[o]=best; }
      if(A[0]>=INF&&A[1]>=INF){fprintf(stderr,"uncovered pair\n");exit(3);}
      loc[mn(A[0],A[1])]++;
      if((A[0]<INF&&A[0]>=0xFFFF)||(A[1]<INF&&A[1]>=0xFFFF)){fprintf(stderr,"16-bit overflow\n");exit(4);} /* max finite A_g is 61242 */
      uint32_t a0=A[0]>=INF?0xFFFF:(uint32_t)A[0], a1=A[1]>=INF?0xFFFF:(uint32_t)A[1];
      key[r1*QQ+r0]=(a0<<16)|a1; }
    #pragma omp critical
    for(int c=0;c<CM;c++) h2[c]+=loc[c];
    free(loc); }
  FILE*f2=fopen(argv[2],"w"); fprintf(f2,"# cost\tmultiplicity\n"); for(int c=0;c<CM;c++) if(h2[c]) fprintf(f2,"%d\t%llu\n",c,(unsigned long long)h2[c]); fclose(f2);
  fprintf(stderr,"k2 written; sorting A keys\n");
  qsort(key,Q2,4,cmp);
  long nu=0; for(long i=0;i<Q2;i++) if(i==0||key[i]!=key[i-1]) key[nu++]=key[i]; /* compact in place after counting below */
  /* recount multiplicities: redo with second pass needs counts; do counting during compaction instead */
  free(key);
  /* second, count-preserving pass (recompute keys; cheaper than storing counts during first compaction) */
  key=malloc(Q2*4);
  #pragma omp parallel for schedule(static)
  for(long r1=0;r1<QQ;r1++) for(long r0=0;r0<QQ;r0++){
      long A[2];
      for(int o=0;o<2;o++){ long best=INF; for(int c1=0;c1<2;c1++) best=mn(best,e[0][c1][r0]+e[c1][o][r1]); A[o]=best; }
      if((A[0]<INF&&A[0]>=0xFFFF)||(A[1]<INF&&A[1]>=0xFFFF)){fprintf(stderr,"16-bit overflow\n");exit(4);} /* max finite A_g is 61242 */
      uint32_t a0=A[0]>=INF?0xFFFF:(uint32_t)A[0], a1=A[1]>=INF?0xFFFF:(uint32_t)A[1];
      key[r1*QQ+r0]=(a0<<16)|a1; }
  qsort(key,Q2,4,cmp);
  uint32_t *uk=malloc(nu*4); uint64_t *un=calloc(nu,8); long t=-1;
  for(long i=0;i<Q2;i++){ if(i==0||key[i]!=key[i-1]) uk[++t]=key[i]; un[t]++; }
  free(key);
  if(t+1!=nu){fprintf(stderr,"pair count mismatch\n");exit(4);}
  static uint32_t tk[QQ]; for(int r2=0;r2<QQ;r2++){ long T0=mn(e[0][0][r2],e[0][1][r2]),T1=mn(e[1][0][r2],e[1][1][r2]);
    uint32_t t0=T0>=INF?0xFFFF:T0,t1=T1>=INF?0xFFFF:T1; tk[r2]=(t0<<16)|t1; }
  qsort(tk,QQ,4,cmp); static uint32_t tu[QQ]; static uint64_t tn[QQ]; int nt=0;
  for(int i=0;i<QQ;i++){ if(i==0||tk[i]!=tk[i-1]) tu[nt++]=tk[i]; tn[nt-1]++; }
  fprintf(stderr,"A-pairs=%ld T-pairs=%d\n",nu,nt);
  uint64_t *h3=calloc(CM,8);
  #pragma omp parallel
  { uint64_t *loc=calloc(CM,8);
    #pragma omp for schedule(dynamic,1)
    for(int j=0;j<nt;j++){ long t0=tu[j]>>16,t1=tu[j]&0xFFFF; uint64_t wt=tn[j];
      for(long i=0;i<nu;i++){ long a0=uk[i]>>16,a1=uk[i]&0xFFFF;
        long x0=(a0==0xFFFF||t0==0xFFFF)?INF:a0+t0, x1=(a1==0xFFFF||t1==0xFFFF)?INF:a1+t1; long c=mn(x0,x1);
        if(c>=INF){fprintf(stderr,"uncovered\n");exit(3);} loc[c]+=wt*un[i]; } }
    #pragma omp critical
    for(int c=0;c<CM;c++) h3[c]+=loc[c];
    free(loc); }
  unsigned long long tot=0; FILE*f=fopen(argv[1],"w"); fprintf(f,"# cost\tmultiplicity\n");
  for(int c=0;c<CM;c++) if(h3[c]){ fprintf(f,"%d\t%llu\n",c,(unsigned long long)h3[c]); tot+=h3[c]; }
  fclose(f); fprintf(stderr,"total=%llu q^3=%llu\n",tot,(unsigned long long)Q2*QQ);
  return 0;
}
