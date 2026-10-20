/* k3_pairs_B.c — Builder B for the k=3 min-cost residue histogram (mod q^3).
   c3(r0 + q R) = min( kappa(r0) + c2[R],  kappa(r0-q) + c2[(R+1) mod q^2] ),  r0 in [0,q), R in [0,q^2),
   where c2 is the full k=2 cost array (same rule as publish/c3a-finite/src/finite_k2.c) and kappa = INF outside D.
   Step 1: build c2; step 2: distinct pairs (c2[R], c2[R+1]) with multiplicities (sorted packed keys);
   step 3: distinct pairs (kappa(r0), kappa(r0-q)); step 4: exact combine -> hist3[c] (uint64).
   Usage: k3_pairs_B OUT_HIST.tsv [pairs_only] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
enum { B=26972, QB=27022, INF=1000000000 };
static const int GENS[]={1971,2016,2100,2628,2688,2800};
static int kap[2*B+1];
static int KA(int d){ return (d<-B||d>B)?INF:kap[d+B]; }
static int cmpu32(const void*a,const void*b){ uint32_t x=*(uint32_t*)a,y=*(uint32_t*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  static unsigned char inM[B+1]; inM[0]=1;
  for(int n=0;n<=B;n++) if(inM[n]) for(int i=0;i<6;i++) if(n+GENS[i]<=B) inM[n+GENS[i]]=1;
  int *M=malloc(sizeof(int)*(B+1)),nM=0; for(int n=0;n<=B;n++) if(inM[n]) M[nM++]=n;
  for(int i=0;i<=2*B;i++) kap[i]=INF;
  for(int i=0;i<nM;i++) for(int j=0;j<nM;j++){ int d=M[i]-M[j]+B,c=M[i]+M[j]; if(c<kap[d]) kap[d]=c; }
  int *Dv=malloc(sizeof(int)*(2*B+1)),*Dc=malloc(sizeof(int)*(2*B+1)),nD=0;
  for(int d=-B;d<=B;d++) if(KA(d)<INF){ Dv[nD]=d; Dc[nD]=KA(d); nD++; }
  fprintf(stderr,"|M|=%d |D|=%d\n",nM,nD);
  long Q2=(long)QB*QB;
  uint16_t *c2=malloc(Q2*2); memset(c2,0xFF,Q2*2);
  #pragma omp parallel for schedule(dynamic,64)
  for(int i=0;i<nD;i++){ long base=Dv[i]; if(base<0) base+=Q2;
    for(int j=0;j<nD;j++){ long r=(base+(long)Dv[j]*QB)%Q2; if(r<0) r+=Q2; int c=Dc[i]+Dc[j];
      /* racy min is unsafe; use atomic CAS loop on 16-bit via 32-bit word */
      uint16_t old=c2[r]; while(c<old){ if(__atomic_compare_exchange_n(&c2[r],&old,(uint16_t)c,0,__ATOMIC_RELAXED,__ATOMIC_RELAXED)) break; }
    } }
  long unc=0; for(long r=0;r<Q2;r++) if(c2[r]==0xFFFF) unc++;
  fprintf(stderr,"c2 built, uncovered=%ld c2[0]=%d c2[h]=%d\n",unc,c2[0],c2[Q2/2]);
  uint32_t *key=malloc(Q2*4);
  for(long R=0;R<Q2;R++) key[R]=((uint32_t)c2[R]<<16)|c2[(R+1)%Q2];
  free(c2);
  qsort(key,Q2,4,cmpu32);
  long nB=0; for(long i=0;i<Q2;i++) if(i==0||key[i]!=key[i-1]) nB++;
  uint32_t *bk=malloc(nB*4); uint64_t *bn=calloc(nB,8); long t=-1;
  for(long i=0;i<Q2;i++){ if(i==0||key[i]!=key[i-1]){ t++; bk[t]=key[i]; } bn[t]++; }
  free(key);
  /* a-pairs */
  uint32_t *ak=malloc(QB*4); for(int r0=0;r0<QB;r0++){ int a0=KA(r0),a1=KA(r0-QB); if(a0>=INF)a0=0xFFFF; if(a1>=INF)a1=0xFFFF; ak[r0]=((uint32_t)a0<<16)|a1; }
  qsort(ak,QB,4,cmpu32); int nA=0; uint32_t *ak2=malloc(QB*4); uint64_t *an=calloc(QB,8);
  for(int i=0;i<QB;i++){ if(i==0||ak[i]!=ak[i-1]){ ak2[nA++]=ak[i]; } an[nA-1]++; }
  fprintf(stderr,"distinct B-pairs=%ld distinct A-pairs=%d\n",nB,nA);
  if(argc>2) return 0;
  enum{CMAX=200000};
  int nth=64; uint64_t *H=calloc((size_t)CMAX*nth,8);
  #pragma omp parallel for schedule(dynamic,1024)
  for(long i=0;i<nB;i++){
    int tid=omp_get_thread_num(); uint64_t *h=H+(size_t)tid*CMAX;
    int b0=bk[i]>>16,b1=bk[i]&0xFFFF; uint64_t w=bn[i];
    for(int j=0;j<nA;j++){ int a0=ak2[j]>>16,a1=ak2[j]&0xFFFF;
      int x0=(a0==0xFFFF||b0==0xFFFF)?INF:a0+b0, x1=(a1==0xFFFF||b1==0xFFFF)?INF:a1+b1; int c=x0<x1?x0:x1;
      if(c>=INF){ fprintf(stderr,"uncovered residue class!\n"); exit(3);} h[c]+=w*an[j]; }
  }
  uint64_t tot=0; FILE*f=fopen(argv[1],"w"); fprintf(f,"# cost\tmultiplicity\n");
  for(int c=0;c<CMAX;c++){ uint64_t s=0; for(int k=0;k<nth;k++) s+=H[(size_t)k*CMAX+c]; if(s){ fprintf(f,"%d\t%llu\n",c,(unsigned long long)s); tot+=s; } }
  fclose(f); fprintf(stderr,"total=%llu (q^3=%llu)\n",(unsigned long long)tot,(unsigned long long)Q2*QB);
  return 0;
}
