// Independent second histogram builder (from hostile skeptic #4), O(q^2) memory-light; writes my_hist.tsv + sumset.txt in cwd.
// Credit: Sylla / xSyllaa.
// independent k=2 min-cost histogram, base q=27022, mask M=<gens> ∩ [0,B]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define B 26972
#define Q 27022
int main(){
  int g[6]={1971,2016,2100,2628,2688,2800};
  static char inM[B+1]; inM[0]=1;
  for(int n=1;n<=B;n++) for(int i=0;i<6;i++) if(n>=g[i]&&inM[n-g[i]]) {inM[n]=1;break;}
  int *M=malloc(sizeof(int)*(B+1)); int m=0; for(int n=0;n<=B;n++) if(inM[n]) M[m++]=n;
  static char inS[2*B+1]; for(int i=0;i<m;i++) for(int j=0;j<m;j++) inS[M[i]+M[j]]=1;
  int ns=0; for(int s=0;s<=2*B;s++) ns+=inS[s];
  fprintf(stderr,"|M|=%d max=%d |M+M|=%d\n",m,M[m-1],ns);
  FILE*fs=fopen("sumset.txt","w"); for(int s=0;s<=2*B;s++) if(inS[s]) fprintf(fs,"%d\n",s); fclose(fs);
  // kappa(d) for d in [-B,B]: min a+b with a-b=d
  const int INF=1<<29; int *ka=malloc(sizeof(int)*(2*B+1)); for(int i=0;i<=2*B;i++) ka[i]=INF;
  for(int i=0;i<m;i++) for(int j=0;j<m;j++){int d=M[i]-M[j]+B; int c=M[i]+M[j]; if(c<ka[d]) ka[d]=c;}
  #define KA(d) (((d)<-B||(d)>B)?INF:ka[(d)+B])
  // per-digit residue cost kb(t)= min(KA(t),KA(t-Q)) for t in [0,Q)
  int *kb=malloc(sizeof(int)*Q); int unc=0;
  for(int t=0;t<Q;t++){int a=KA(t),b=KA(t-Q); kb[t]=a<b?a:b; if(kb[t]>=INF) unc++;}
  fprintf(stderr,"uncovered digit residues=%d kappa(0)=%d\n",unc,KA(0));
  // block residue r=r0+r1*Q (mod Q^2): lifts e=d0+d1*Q, d0 in {r0,r0-Q}, d1 ≡ r1 (+1 if d0=r0-Q)
  long long *h=calloc(2*2*B+2,sizeof(long long)); long long tot=0; int mx=0; long long unc2=0;
  int *c2row=malloc(sizeof(int)*Q);
  long long asym=0; int ch=-1;
  for(int r1=0;r1<Q;r1++){
    for(int r0=0;r0<Q;r0++){
      int a=KA(r0); a = (a>=INF)?INF: a + kb[r1];
      int b=KA(r0-Q); b = (b>=INF)?INF: b + kb[(r1+1)%Q];
      int c=a<b?a:b;
      if(c>=INF){unc2++;continue;}
      h[c]++; tot++; if(c>mx) mx=c;
      if(r0==0 && r1==Q/2) ch=c;
    }
  }
  // symmetry check: c(-r) vs c(r) on sample
  srand(1); int bad=0;
  for(int s=0;s<200000;s++){ long long r=((long long)rand()*7919+rand())%((long long)Q*Q); long long nr=((long long)Q*Q-r)%((long long)Q*Q);
    int cc[2]; long long rr[2]={r,nr};
    for(int k=0;k<2;k++){int r0=rr[k]%Q,r1=rr[k]/Q; int a=KA(r0); a=(a>=INF)?INF:a+kb[r1]; int b=KA(r0-Q); b=(b>=INF)?INF:b+kb[(r1+1)%Q]; cc[k]=a<b?a:b;}
    if(cc[0]!=cc[1]) bad++; }
  fprintf(stderr,"tot=%lld uncovered=%lld max=%d c(h)=%d symmetry_bad=%d\n",tot,unc2,mx,ch,bad);
  FILE*f=fopen("my_hist.tsv","w"); fprintf(f,"# cost\tmultiplicity\n"); for(int c=0;c<=mx;c++) if(h[c]) fprintf(f,"%d\t%lld\n",c,h[c]); fclose(f);
  return 0;
}
