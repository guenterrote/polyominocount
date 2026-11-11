#define MAXL 160
#define PSELF 0xFF
#define PGHOST 0xFE
#define MAXN 2*(MAXL+4) 
#define GH(i) (GOFF+(i) ) 
#define ISGHOST(a) ((a) > GOFF) 
/*2:*/
#line 104 "diag180.w"

/*3:*/
#line 111 "diag180.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <stdint.h> 
#include <time.h> 

/*:3*/
#line 105 "diag180.w"

/*14:*/
#line 350 "diag180.w"

typedef struct{
int f[MAXL];
int partner[MAXL];
int parent[MAXN];
int next;
}work_t;

/*:14*//*22:*/
#line 483 "diag180.w"

enum{OK,DEAD,COMPLETE};

/*:22*//*27:*/
#line 597 "diag180.w"

enum{C_EMPTY,C_CONTINUE,C_SINGLETON,C_FIRST,C_MIDDLE,C_LAST};

/*:27*//*32:*/
#line 718 "diag180.w"

#ifdef CNT128
typedef unsigned __int128 cnt_t;
#define CNT_BITS "128"
#else
#line 723 "diag180.w"
 typedef uint64_t cnt_t;
#define CNT_BITS "64"
#endif
#line 726 "diag180.w"

/*:32*//*36:*/
#line 769 "diag180.w"

typedef struct{
cnt_t c;
uint32_t next;
uint8_t wt;
}rec_t;

typedef struct{
size_t n,cap,hmask;
uint64_t*keys;
uint32_t*head;
rec_t*rec;
size_t nrec,rcap;
uint32_t*hash;
}table_t;

/*:36*/
#line 106 "diag180.w"

/*8:*/
#line 210 "diag180.w"

static int NMAX;
static int L,vlo,T,S,E;
static int NW;

/*:8*//*13:*/
#line 337 "diag180.w"

static int MAXC;
static int KL;

/*:13*//*15:*/
#line 358 "diag180.w"

static int GOFF;

/*:15*//*26:*/
#line 594 "diag180.w"

static int KW;

/*:26*//*33:*/
#line 727 "diag180.w"

static uint64_t MOD= 0;

/*:33*//*53:*/
#line 1173 "diag180.w"

static int CUR_U,USTAR,U0;
static int HEIGHT= 1;
static int H0,S0;
static int NOPT,OPTD[2*MAXL+4],OPTP[2*MAXL+4],OPTS[2*MAXL+4];
static int NCL,CLD[MAXL],CLP[MAXL],CLS[MAXL];
static int DEDGE;
static int CUR_TCH;

/*:53*//*69:*/
#line 1574 "diag180.w"

static cnt_t*GF;
static cnt_t*LT,*EQ;
static size_t MAXSTATES,MAXBYTES,MAXRECS,MAXKW;

/*:69*//*71:*/
#line 1599 "diag180.w"

static int VERB= 1;
static double PSEC= 60;
static int G_C,G_D,G_VHI;
static double T_START,T_LAST;

/*:71*//*83:*/
#line 1901 "diag180.w"

static int REFLECT= 1;
static int OLDBOUND= 0;
static const char*RESFILE= "results180.py";
static int LOOKAHEAD= 1;
static int CHECK= 0;

/*:83*/
#line 107 "diag180.w"

/*9:*/
#line 233 "diag180.w"

static int uspan(int k){return 2*k+1-E;}

/*:9*//*16:*/
#line 362 "diag180.w"

static inline int find(work_t*w,int a)
{
while(w->parent[a]!=a){
w->parent[a]= w->parent[w->parent[a]];
a= w->parent[a];
}
return a;
}

/*:16*//*17:*/
#line 374 "diag180.w"

static inline int rho(work_t*w,int a)
{
a= find(w,a);
if(ISGHOST(a))
return find(w,a-GOFF);
int p= w->partner[a];
if(p==-1)
return a;
if(p==-2)
return find(w,GH(a));
return find(w,p);
}

/*:17*//*18:*/
#line 391 "diag180.w"

static inline void link_(work_t*w,int a,int b)
{
a= find(w,a);
b= find(w,b);
if(a==b)
return;
if(ISGHOST(a)){
int x= a;
a= b;
b= x;
}
w->parent[b]= a;
}

/*:18*//*19:*/
#line 410 "diag180.w"

static inline void unite(work_t*w,int a,int b)
{
int ra= rho(w,a),rb= rho(w,b);
link_(w,a,b);
link_(w,ra,rb);
}

/*:19*//*20:*/
#line 419 "diag180.w"

static void decode(work_t*w,const uint8_t*key)
{
int m= 0;
for(int i= 0;i<L;i++){
w->f[i]= key[i];
if(key[i]> m)
m= key[i];
}
for(int a= 0;a<MAXN;a++)
w->parent[a]= a;
for(int j= 1;j<=MAXC;j++){
int q= key[L+j-1];
w->partner[j]= q==PSELF?-1:q==PGHOST?-2:q;
if(q&&j> m)
m= j;
}
w->next= m+1;
}

/*:20*//*21:*/
#line 442 "diag180.w"

static int encode(work_t*w,uint8_t*key,int tch)
{
int map[MAXN];
int rep[MAXL];
int m= 0;
memset(map,0,sizeof(map));
memset(key,0,KL);
for(int i= 0;i<L;i++){
if(!w->f[i])
continue;
int r= find(w,w->f[i]);
if(!map[r]){
map[r]= ++m;
rep[m]= r;
}
key[i]= map[r];
}
for(int j= 1;j<=m;j++){
int q= rho(w,rep[j]);
key[L+j-1]= q==rep[j]?PSELF:map[q]?map[q]:PGHOST;
}
key[KL-1]= tch;
return m;
}

/*:21*//*23:*/
#line 486 "diag180.w"

static int retire(work_t*w,int i)
{
int lab= w->f[i];
if(!lab)
return OK;
w->f[i]= 0;
int r= find(w,lab);
int any= 0;
for(int j= 0;j<L;j++)
if(w->f[j]){
any= 1;
if(find(w,w->f[j])==r)
return OK;
}
int p= rho(w,r);
if(p==r)
return any?DEAD:COMPLETE;
if(ISGHOST(p))
return DEAD;
w->partner[find(w,p)]= -2;
return OK;
}

/*:23*//*24:*/
#line 526 "diag180.w"

static int step(work_t*w,const uint8_t*sk,int tch,int i,int occ,int last,
int first_twin,uint8_t*key,int*ntp,int*mp)
{
decode(w,sk);
if(occ){
int nw= w->next++;
w->partner[nw]= -2;
w->f[i]= nw;
for(int nb= i-1;nb<=i+1;nb+= 2)
if(nb>=0&&nb<L&&w->f[nb])
unite(w,nw,w->f[nb]);
if(first_twin){
for(int jj= L-2-i;jj<=L-i;jj+= 2){
if(jj==i)
unite(w,nw,rho(w,nw));
else if(jj>=0&&jj<L&&jj<i&&w->f[jj])
unite(w,nw,rho(w,w->f[jj]));
}
}
}
int status= OK;
if(i-1>=0)
status= retire(w,i-1);
if(status==OK&&last&&i+1<L)
status= retire(w,i+1);
if(status!=OK)
return status;
*ntp= tch||(occ&&(i==0||i==L-1));
*mp= encode(w,key,*ntp);
return OK;
}

/*:24*//*29:*/
#line 621 "diag180.w"

static inline void putbits(uint64_t*w,int*pos,uint64_t v,int nb)
{
int i= *pos>>6,o= *pos&63;
w[i]|= v<<o;
if(o+nb> 64)
w[i+1]|= v>>(64-o);
*pos+= nb;
}

static inline int getbits(const uint64_t*w,int*pos,int nb)
{
int i= *pos>>6,o= *pos&63;
uint64_t v= w[i]>>o;
if(o+nb> 64)
v|= w[i+1]<<(64-o);
*pos+= nb;
return v&((1u<<nb)-1);
}

/*:29*//*30:*/
#line 643 "diag180.w"

static void pack(const uint8_t*key,uint64_t*w)
{
int first[MAXL+2],lastrun[MAXL+2];
int pos= 0,m= 0;
memset(w,0,KW*8);
for(int i= 0;i<L;i++)
if(key[i]&&(i==0||!key[i-1])){
int l= key[i];
if(l> m){
m= l;
first[l]= i;
}
lastrun[l]= i;
}
for(int i= 0;i<L;i++){
int l= key[i],code;
if(!l)
code= C_EMPTY;
else if(i> 0&&key[i-1])
code= C_CONTINUE;
else if(first[l]==i)
code= lastrun[l]==i?C_SINGLETON:C_FIRST;
else
code= lastrun[l]==i?C_LAST:C_MIDDLE;
putbits(w,&pos,code,3);
}
for(int j= 1;j<=m;j++)
putbits(w,&pos,key[L+j-1]!=PGHOST,1);
putbits(w,&pos,key[KL-1],1);
}

/*:30*//*31:*/
#line 677 "diag180.w"

static void unpack(const uint64_t*w,uint8_t*key)
{
int stack[MAXL+2],sp= 0,m= 0,pos= 0;
memset(key,0,KL);
for(int i= 0;i<L;i++){
int code= getbits(w,&pos,3);
switch(code){
case C_EMPTY:key[i]= 0;break;
case C_CONTINUE:key[i]= key[i-1];break;
case C_SINGLETON:key[i]= ++m;break;
case C_FIRST:key[i]= ++m;stack[sp++]= m;break;
case C_MIDDLE:key[i]= stack[sp-1];break;
case C_LAST:key[i]= stack[--sp];break;
}
}
sp= 0;
for(int j= 1;j<=m;j++)
if(getbits(w,&pos,1))
stack[sp++]= j;
else
key[L+j-1]= PGHOST;
for(int a= 0,b= sp-1;a<=b;a++,b--){
key[L+stack[a]-1]= a==b?PSELF:stack[b];
key[L+stack[b]-1]= a==b?PSELF:stack[a];
}
key[KL-1]= getbits(w,&pos,1);
}

/*:31*//*34:*/
#line 730 "diag180.w"

static inline cnt_t addm(cnt_t a,cnt_t b)
{
if(!MOD)
return a+b;
cnt_t s= a+b;
if(s<a||s>=MOD)
s-= MOD;
return s;
}

/*:34*//*35:*/
#line 745 "diag180.w"

static const char*cs(cnt_t x)
{
static char buf[8][44];
static int r;
char*e= buf[r= (r+1)&7]+43,*p= e;
*p= 0;
do
*--p= '0'+(int)(x%10),x/= 10;
while(x);
return p;
}

/*:35*//*37:*/
#line 785 "diag180.w"

static void*xrealloc(void*p,size_t s)
{
p= realloc(p,s);
if(!p){
fprintf(stderr,"out of memory\n");
exit(1);
}
return p;
}

/*:37*//*38:*/
#line 796 "diag180.w"

static void tab_init(table_t*t,size_t cap)
{
t->n= 0;
t->cap= cap;
size_t h= 16;
while(h<2*cap)
h<<= 1;
t->hmask= h-1;
t->keys= xrealloc(0,cap*KW*8);
t->head= xrealloc(0,cap*4);
t->hash= calloc(h,sizeof(uint32_t));
t->rcap= 2*cap;
t->nrec= 0;
t->rec= xrealloc(0,t->rcap*sizeof(rec_t));
}

static void tab_free(table_t*t)
{
free(t->keys);
free(t->head);
free(t->rec);
free(t->hash);
}

static size_t tab_bytes(const table_t*t)
{
return t->n*(KW*8+4)+t->nrec*sizeof(rec_t)+(t->hmask+1)*4;
}

/*:38*//*39:*/
#line 827 "diag180.w"

static inline uint64_t hkey(const uint64_t*k)
{
uint64_t h= 0x9E3779B97F4A7C15ull;
for(int i= 0;i<KW;i++){
h^= k[i];
h*= 0xBF58476D1CE4E5B9ull;
h^= h>>31;
}
return h;
}

/*:39*//*40:*/
#line 841 "diag180.w"

static void tab_grow(table_t*t)
{
size_t cap= t->cap*2;
t->keys= xrealloc(t->keys,cap*KW*8);
t->head= xrealloc(t->head,cap*4);
size_t h= (t->hmask+1)*2;
free(t->hash);
t->hash= calloc(h,sizeof(uint32_t));
if(!t->hash){
fprintf(stderr,"out of memory\n");
exit(1);
}
t->hmask= h-1;
t->cap= cap;
for(size_t i= 0;i<t->n;i++){
size_t p= hkey(t->keys+i*KW)&t->hmask;
while(t->hash[p])
p= (p+1)&t->hmask;
t->hash[p]= i+1;
}
}

static void tab_clear(table_t*t)
{
t->n= 0;
t->nrec= 0;
memset(t->hash,0,(t->hmask+1)*sizeof(uint32_t));
}

/*:40*//*41:*/
#line 874 "diag180.w"

static size_t tab_get(table_t*t,const uint8_t*key)
{
uint64_t w[12];
pack(key,w);
if(CHECK){
uint8_t k2[2*MAXL+4];
unpack(w,k2);
if(memcmp(k2,key,KL)){
fprintf(stderr,"pack/unpack mismatch\n");
exit(2);
}
}
size_t p= hkey(w)&t->hmask;
while(t->hash[p]){
size_t i= t->hash[p]-1;
if(!memcmp(t->keys+i*KW,w,KW*8))
return i;
p= (p+1)&t->hmask;
}
if(t->n==t->cap){
tab_grow(t);
p= hkey(w)&t->hmask;
while(t->hash[p])
p= (p+1)&t->hmask;
}
size_t i= t->n++;
t->hash[p]= i+1;
memcpy(t->keys+i*KW,w,KW*8);
t->head[i]= 0;
return i;
}

/*:41*//*42:*/
#line 908 "diag180.w"

static void tab_add(table_t*t,size_t si,int wt,cnt_t c)
{
for(uint32_t r= t->head[si];r;r= t->rec[r-1].next)
if(t->rec[r-1].wt==wt){
t->rec[r-1].c= addm(t->rec[r-1].c,c);
return;
}
if(t->nrec==t->rcap){
t->rcap*= 2;
t->rec= xrealloc(t->rec,t->rcap*sizeof(rec_t));
}
rec_t*r= &t->rec[t->nrec++];
r->c= c;
r->wt= wt;
r->next= t->head[si];
t->head[si]= t->nrec;
}

/*:42*//*43:*/
#line 930 "diag180.w"

static void tab_load(const table_t*t,size_t si,cnt_t*gf)
{
for(uint32_t r= t->head[si];r;r= t->rec[r-1].next)
gf[t->rec[r-1].wt]= t->rec[r-1].c;
}

static void tab_unload(const table_t*t,size_t si,cnt_t*gf)
{
for(uint32_t r= t->head[si];r;r= t->rec[r-1].next)
gf[t->rec[r-1].wt]= 0;
}

/*:43*//*44:*/
#line 949 "diag180.w"

static void insert_pruned(table_t*t,const uint8_t*key,const cnt_t*gf,
int add,int lb)
{
int lim= NMAX-lb-add;
int any= 0;
for(int wt= 0;wt<=lim;wt++)
if(gf[wt]){
any= 1;
break;
}
if(!any)
return;
size_t si= tab_get(t,key);
for(int wt= 0;wt<=lim;wt++)
if(gf[wt])
tab_add(t,si,wt+add,gf[wt]);
}

/*:44*//*46:*/
#line 982 "diag180.w"

static int completion_possible(const uint8_t*key)
{
int m= 0;
for(int j= 0;j<MAXC;j++){
int q= key[L+j];
if(!q)
break;
m++;
if(q!=PGHOST)
return 1;
}
return m==0;
}

/*:46*//*54:*/
#line 1184 "diag180.w"

static inline int frow(int j)
{
return((vlo+j-CUR_U)%2+2)%2?CUR_U-1:CUR_U;
}

/*:54*//*55:*/
#line 1190 "diag180.w"

static inline void reach(int h)
{
if(h> H0)
H0= h;
}

static inline void sreach(int sc)
{
if(sc<S0)
S0= sc;
}

/*:55*//*56:*/
#line 1205 "diag180.w"

static inline int both_cost(int h,int j)
{
int e= j<L-1-j?j:L-1-j;
int r= USTAR-h;
return e> r?e:r;
}

/*:56*//*57:*/
#line 1216 "diag180.w"

static int gap_score(int a,int ra,int e,int re_)
{
int best= 1<<20;
for(int j= a+1;j<=e;j++){
int t1= ra+1+(j-a-1),t2= re_+(e-j);
int c= both_cost(t1<t2?t1:t2,j);
if(c<best)
best= c;
}
return best;
}

/*:57*//*58:*/
#line 1235 "diag180.w"

static int cycle(const int*g,const int*pk,const int*sc,int n,int*m1p,int*p1p,int*s1p)
{
int sum= 0,m1= -1,n1= 0,m2= 0,p1= -1000,s1= 1<<20;
for(int t= 0;t<n;t++){
sum+= g[t];
if(g[t]> m1){
if(m1>=0)
m2= m1;
m1= g[t];
n1= 1;
p1= pk[t];
s1= sc[t];
}else if(g[t]==m1){
n1++;
if(pk[t]> p1)
p1= pk[t];
if(sc[t]<s1)
s1= sc[t];
}else if(g[t]> m2)
m2= g[t];
}
for(int t= 0;t<n;t++)
if(n1>=2||g[t]!=m1){
reach(pk[t]);
sreach(sc[t]);
}
if(n1==1){
OPTD[NOPT]= m1-m2;
OPTS[NOPT]= s1;
OPTP[NOPT++]= p1;
}
if(m1p){
*m1p= m1;
*p1p= p1;
*s1p= s1;
}
return sum-m1;
}

/*:58*//*59:*/
#line 1278 "diag180.w"

static inline int finish(int h,int s)
{
int top= USTAR-h> 0?USTAR-h:0;
if(CUR_TCH)
return HEIGHT?top:0;
if(!HEIGHT)
return DEDGE;
int sep= DEDGE+top;
return s<sep?s:sep;
}

/*:59*//*60:*/
#line 1292 "diag180.w"

static int height_extra(int h,int s)
{
int ex= finish(h,s);
for(int o= 0;o<NOPT&&ex> 0;o++){
int e2= OPTD[o]+finish(OPTP[o]> h?OPTP[o]:h,OPTS[o]<s?OPTS[o]:s);
if(e2<ex)
ex= e2;
}
return ex;
}

/*:60*//*61:*/
#line 1307 "diag180.w"

static int cactus_bound(const uint8_t*key,int tch,int kink)
{
/*62:*/
#line 1321 "diag180.w"

int rs[MAXL],re[MAXL],rl[MAXL],R= 0;
int lastrun[MAXL+2];
for(int i= 0;i<L;i++)
if(key[i]){
if(i==0||!key[i-1]){
rs[R]= i;
rl[R]= key[i];
R++;
}
re[R-1]= i;
lastrun[key[i]]= R-1;
}
if(!R)
return 0;

/*:62*/
#line 1310 "diag180.w"

/*63:*/
#line 1337 "diag180.w"

H0= -1000;
S0= 1<<20;
NOPT= 0;
NCL= 0;
CUR_TCH= tch;
for(int i= 0;i<L;i++)
if(key[i]){
reach(frow(i));
sreach(both_cost(frow(i),i));
}

/*:63*/
#line 1311 "diag180.w"

int cost= 0;
/*64:*/
#line 1356 "diag180.w"


int sl[MAXL],sb[MAXL],sp= 0;
int gb[MAXL],pb[MAXL],qb[MAXL],nb= 0;
int toplab[MAXL],topgap[MAXL],toppk[MAXL],topsc[MAXL],nt= 0;
for(int j= 0;j<R;j++){
int X= rl[j];
if(j> 0){
int a= re[j-1],b= rs[j];
int g= b-a-1+(b==kink);

int pk= (g-1+(b==kink)+frow(a)+frow(b)+2)/2;
int sc= tch?1<<20:gap_score(a,frow(a),b==kink?b+1:b-1,frow(b)+1);
if(sp){
gb[nb]= g;
qb[nb]= sc;
pb[nb++]= pk;
}else{
topgap[nt-1]= g;
topsc[nt-1]= sc;
toppk[nt-1]= pk;
}
}
if(sp&&sl[sp-1]==X){
sp--;
cost+= cycle(gb+sb[sp],pb+sb[sp],qb+sb[sp],nb-sb[sp],0,0,0);
nb= sb[sp];
}else if(!sp)
toplab[nt++]= X;
if(j<lastrun[X]){
sl[sp]= X;
sb[sp]= nb;
sp++;
}
}

/*:64*/
#line 1313 "diag180.w"

/*65:*/
#line 1395 "diag180.w"


int cr[MAXL],m= 0;
for(int t= 0;t<nt;t++)
if(key[L+toplab[t]-1]!=PGHOST)
cr[m++]= t;
if(m){
for(int t= 0;t<cr[0];t++){
cost+= topgap[t];
reach(toppk[t]);
sreach(topsc[t]);
}
for(int t= cr[m-1];t<nt-1;t++){
cost+= topgap[t];
reach(toppk[t]);
sreach(topsc[t]);
}
int k= m/2;
NCL= 0;
for(int i= 0;i<k;i++){
int g[2*MAXL],pk[2*MAXL],scs[2*MAXL],n= 0,m1,p1,s1;
int a= cr[i],b= (i<k-1)?cr[i+1]:cr[m-1-i];
for(int t= a;t<b;t++){
g[n]= topgap[t];
scs[n]= topsc[t];
pk[n++]= toppk[t];
}
if(i<k-1)
for(int t= cr[m-2-i];t<cr[m-1-i];t++){
g[n]= topgap[t];
scs[n]= topsc[t];
pk[n++]= toppk[t];
}
cost+= cycle(g,pk,scs,n,&m1,&p1,&s1);
if(i==k-1){



CLD[NCL]= m1;
CLS[NCL]= s1;
CLP[NCL++]= p1;
}
}
if(m&1)
NCL= 0;
}else{
for(int t= 0;t<nt-1;t++){
cost+= topgap[t];
reach(toppk[t]);
sreach(topsc[t]);
}
}

/*:65*/
#line 1314 "diag180.w"

/*66:*/
#line 1448 "diag180.w"

DEDGE= 0;
if(!tch){
int d1= rs[0],d2= L-1-re[R-1];
DEDGE= d1<d2?d1:d2;
}


if(NCL){
int ex= 1<<20;
for(int c= 0;c<NCL;c++){
int e2= CLD[c]+height_extra(CLP[c]> H0?CLP[c]:H0,CLS[c]<S0?CLS[c]:S0);
if(e2<ex)
ex= e2;
}
cost+= ex;
}else
cost+= height_extra(H0,S0);

/*:66*/
#line 1315 "diag180.w"

return cost;
}

/*:61*//*67:*/
#line 1471 "diag180.w"

static int exact_bound(const uint8_t*key,int tch)
{
int sup[MAXL+1];
int occ[MAXL],no= 0;
for(int j= 1;j<=MAXC;j++){
int q= key[L+j-1];
if(!q)
break;
sup[j]= (q!=PSELF&&q!=PGHOST&&q<j)?q:j;
}
for(int i= 0;i<L;i++)
if(key[i])
occ[no++]= i;
if(!no)
return 0;

int par[MAXL+1];
int ncls= 0;
int seen[MAXL+1];
memset(seen,0,sizeof(seen));
for(int a= 0;a<no;a++){
int s= sup[key[occ[a]]];
if(!seen[s]){
seen[s]= 1;
ncls++;
}
par[s]= s;
}
int cost= 0;
if(ncls> 1){

int ne= no-1;
int wgt[MAXL],ord[MAXL],bucket[MAXL+1];
memset(bucket,0,sizeof(bucket));
for(int a= 0;a<ne;a++){
wgt[a]= occ[a+1]-occ[a]-1;
bucket[wgt[a]]++;
}
for(int g= 1;g<=L;g++)
bucket[g]+= bucket[g-1];
for(int a= ne-1;a>=0;a--)
ord[--bucket[wgt[a]]]= a;
for(int e= 0;e<ne&&ncls> 1;e++){
int a= ord[e];
int x= sup[key[occ[a]]],y= sup[key[occ[a+1]]];
while(par[x]!=x)
x= par[x];
while(par[y]!=y)
y= par[y];
if(x!=y){
par[x]= y;
cost+= wgt[a];
ncls--;
}
}
}
if(!tch){
int d1= occ[0],d2= L-1-occ[no-1];
cost+= d1<d2?d1:d2;
}
return cost;
}

/*:67*//*68:*/
#line 1542 "diag180.w"

static int bound(const uint8_t*key,int nt,int k,int kmin,int i,int last,int first_twin)
{
int lb= kmin> k?2*(kmin-k):0;
if(!nt){
int dist= L;
for(int j= 0;j<L;j++)
if(key[j]){
int dd= j<L-1-j?j:L-1-j;
if(dd<dist)
dist= dd;
}
if(2*dist> lb)
lb= 2*dist;
}
if(!first_twin||last){
if(!completion_possible(key))
return-1;
int kink= (!last&&!key[i]&&i+1<L&&key[i+1])?i+1:-1;
CUR_U= U0+k-1;
USTAR= U0+kmin-1;
int eb= 2*(OLDBOUND?exact_bound(key,nt):cactus_bound(key,nt,kink));
if(eb> lb)
lb= eb;
}
return lb;
}

/*:68*//*70:*/
#line 1581 "diag180.w"

static void record(int ktop,const cnt_t*gf,int tch)
{
int h= uspan(ktop);
if(!tch||h<L)
return;
cnt_t*tgt= h==L?EQ:LT;
for(int wt= 0;wt<=NMAX;wt++)
if(gf[wt])
tgt[wt]= addm(tgt[wt],gf[wt]);
}

/*:70*//*72:*/
#line 1605 "diag180.w"

static double wall(void)
{
struct timespec ts;
clock_gettime(CLOCK_MONOTONIC,&ts);
return ts.tv_sec+1e-9*ts.tv_nsec;
}
static void stat_line(const char*tag,int k,int u,int ci,int nu,
const table_t*cur,const table_t*nxt)
{
double now= wall();
T_LAST= now;
cnt_t done= 0;
for(int n= 0;n<=NMAX;n++)
done= addm(done,addm(LT[n],EQ[n]));
fprintf(stderr,"%s c=%d d=%d vhi=%d L=%d row k=%d u=%d cell %d/%d: "
"%zu states, %zu records, %.1f MB, done %s, %.1fs wall, %.1fs cpu\n",
tag,G_C,G_D,G_VHI,L,k,u,ci,nu,cur->n,cur->nrec,
(tab_bytes(cur)+tab_bytes(nxt))/1e6,cs(done),
now-T_START,(double)clock()/CLOCKS_PER_SEC);
fflush(stderr);
}

/*:72*//*73:*/
#line 1629 "diag180.w"

static void count_case(int c,int d,int vhi)
{
/*74:*/
#line 1650 "diag180.w"

G_C= c,G_D= d,G_VHI= vhi;
S= c+d;
T= c-d;
vlo= T-vhi;
L= vhi-vlo+1;
E= S&1;
memset(LT,0,NW*sizeof(cnt_t));
memset(EQ,0,NW*sizeof(cnt_t));
if(L==1){
if(c==0&&d==0)
EQ[1]= addm(EQ[1],1);
return;
}
if(L> MAXL){
fprintf(stderr,"strip too wide\n");
exit(1);
}
MAXC= (L+1)/2+1;
GOFF= MAXC+1;
KL= L+MAXC+1;
KW= (3*L+MAXC+1+63)/64;
NW= NMAX+1;

int kmin= 0;
while(uspan(kmin)<L)
kmin++;

/*:74*/
#line 1632 "diag180.w"

/*75:*/
#line 1680 "diag180.w"

table_t A,B,*cur= &A,*nxt= &B;
tab_init(&A,1024);
tab_init(&B,1024);
uint8_t key[MAXL*2+4];
cnt_t one[MAXL*2+4];
work_t w;
int u0;

/*:75*/
#line 1633 "diag180.w"

/*76:*/
#line 1696 "diag180.w"

if(E==0){
int um= S/2;
int half[MAXL],nh= 0;
for(int i= 0;i<L;i++)
if(((vlo+i-um)%2+2)%2==0&&2*(vlo+i)>=T)
half[nh++]= i;

for(uint64_t mask= 1;mask<(1ull<<nh);mask++){
int cells[MAXL],nc= 0,isc[MAXL];
memset(isc,0,sizeof(int)*L);
for(int b= 0;b<nh;b++)
if(mask>>b&1){
isc[half[b]]= 1;
isc[L-1-half[b]]= 1;
}
for(int i= 0;i<L;i++)
if(isc[i])
cells[nc++]= i;
if(nc> NMAX)
continue;
int lab[MAXL];
for(int a= 0;a<MAXN;a++)
w.parent[a]= a;
for(int i= 0;i<L;i++)
w.f[i]= 0;
for(int a= 0;a<nc;a++){
lab[cells[a]]= a+1;
w.f[cells[a]]= a+1;
}
for(int a= 0;a<nc;a++){
int i= cells[a],j= L-1-i;
w.partner[a+1]= j==i?-1:lab[j];
}
int tch= isc[0]||isc[L-1];
encode(&w,key,tch);
memset(one,0,sizeof(cnt_t)*NW);
one[nc]= 1;
insert_pruned(cur,key,one,0,0);
}
u0= um+1;
U0= u0;
}else{
memset(key,0,KL);
memset(one,0,sizeof(cnt_t)*NW);
one[0]= 1;
insert_pruned(cur,key,one,0,0);
u0= (S+1)/2;
U0= u0;
}

/*:76*/
#line 1634 "diag180.w"

int cells_u[MAXL];
for(int k= 1;cur->n;k++){
/*77:*/
#line 1749 "diag180.w"

int u= u0+k-1,nu= 0;
for(int i= 0;i<L;i++)
if(((vlo+i-u)%2+2)%2==0)
cells_u[nu++]= i;
int first_twin= (E==1&&k==1);

/*:77*/
#line 1637 "diag180.w"

for(int ci= 0;ci<nu;ci++){
/*78:*/
#line 1756 "diag180.w"

int i= cells_u[ci];
int last= ci==nu-1;
tab_clear(nxt);
for(size_t si= 0;si<cur->n;si++){
/*79:*/
#line 1767 "diag180.w"

uint8_t sk[MAXL*2+4];
unpack(cur->keys+si*KW,sk);
const cnt_t*gf= GF;
tab_load(cur,si,GF);
int tch= sk[KL-1];
int empty= 1;
for(int j= 0;j<L;j++)
if(sk[j]){
empty= 0;
break;
}
for(int occ= 0;occ<=1;occ++){
if(occ&&k> 1&&empty)
continue;
int nt,m;
int status= step(&w,sk,tch,i,occ,last,first_twin,key,&nt,&m);
if(status==DEAD)
continue;
if(status==COMPLETE){
record(k-1,gf,tch);
continue;
}
int add= occ?2:0;
if(m==0){
if(first_twin&&!last)
insert_pruned(nxt,key,gf,add,0);
continue;
}
int lb= bound(key,nt,k,kmin,i,last,first_twin);
if(lb<0)
continue;
/*80:*/
#line 1814 "diag180.w"

if(LOOKAHEAD&&!last&&!first_twin&&(LOOKAHEAD==2||(!key[i]&&i+1<L&&key[i+1]))){

int i2= cells_u[ci+1],last2= ci+1==nu-1;
int la= 1<<20;
uint8_t k2[MAXL*2+4];
for(int o2= 0;o2<=1&&la;o2++){
int nt2,m2;
int st2= step(&w,key,nt,i2,o2,last2,0,k2,&nt2,&m2);
if(st2==DEAD)
continue;

if(st2==COMPLETE&&!(nt&&uspan(k-1)>=L))
continue;
int val= 2*o2;
if(st2==OK){
if(m2==0)
continue;
int lb2= bound(k2,nt2,k,kmin,i2,last2,0);
if(lb2<0)
continue;
val+= lb2;
}
if(val<la)
la= val;
}
if(la> lb)
lb= la;
}

/*:80*/
#line 1799 "diag180.w"

insert_pruned(nxt,key,gf,add,lb);
}
tab_unload(cur,si,GF);

/*:79*/
#line 1761 "diag180.w"

}
/*81:*/
#line 1844 "diag180.w"

table_t*x= cur;
cur= nxt;
nxt= x;
if(cur->n> MAXSTATES)
MAXSTATES= cur->n;
size_t by= tab_bytes(cur)+tab_bytes(nxt);
if(by> MAXBYTES)
MAXBYTES= by,MAXRECS= cur->nrec,MAXKW= KW;
if(!last&&(VERB>=2||(VERB==1&&wall()-T_LAST>=PSEC)))
stat_line("#cell",k,u,ci+1,nu,cur,nxt);

/*:81*/
#line 1763 "diag180.w"


/*:78*/
#line 1639 "diag180.w"

}
/*82:*/
#line 1864 "diag180.w"

if(REFLECT&&T==0&&cur->n){

tab_clear(nxt);
uint8_t rkey[MAXL*2+4];
for(size_t si= 0;si<cur->n;si++){
uint8_t sk[MAXL*2+4];
unpack(cur->keys+si*KW,sk);
decode(&w,sk);
for(int a= 0,b= L-1;a<b;a++,b--){
int x= w.f[a];
w.f[a]= w.f[b];
w.f[b]= x;
}
encode(&w,rkey,sk[KL-1]);
const uint8_t*use= memcmp(rkey,sk,KL)<0?rkey:sk;
tab_load(cur,si,GF);
insert_pruned(nxt,use,GF,0,0);
tab_unload(cur,si,GF);
}
table_t*x= cur;
cur= nxt;
nxt= x;
}

/*:82*/
#line 1641 "diag180.w"

if(VERB>=1)
stat_line("#row ",k,u,nu,nu,cur,nxt);
}
tab_free(&A);
tab_free(&B);
}

/*:73*/
#line 108 "diag180.w"

/*84:*/
#line 1908 "diag180.w"

int main(int argc,char**argv)
{
/*85:*/
#line 1926 "diag180.w"

int a= 1;
while(a<argc&&argv[a][0]=='-'){
if(!strcmp(argv[a],"-m")&&a+1<argc)
MOD= strtoull(argv[++a],0,10),a++;
else if(!strcmp(argv[a],"-o")&&a+1<argc)
RESFILE= argv[++a],a++;
else if(!strcmp(argv[a],"-A"))
LOOKAHEAD= 2,a++;
else if(!strcmp(argv[a],"-N"))
LOOKAHEAD= 0,a++;
else if(!strcmp(argv[a],"-H"))
HEIGHT= 0,a++;
else if(!strcmp(argv[a],"-K"))
OLDBOUND= 1,a++;
else if(!strcmp(argv[a],"-C"))
CHECK= 1,a++;
else if(!strcmp(argv[a],"-v")&&a+1<argc)
VERB= atoi(argv[++a]),a++;
else if(!strcmp(argv[a],"-p")&&a+1<argc)
PSEC= atof(argv[++a]),a++;
else if(!strcmp(argv[a],"-R"))
REFLECT= 0,a++;
else
break;
}
if(a<argc&&argv[a][0]=='-'){
fprintf(stderr,"%s: unknown option %s\n",argv[0],argv[a]);
return 1;
}
if(a>=argc){
fprintf(stderr,"usage: %s [-m M] [-R] N [c d vhi]\n",argv[0]);
return 1;
}
NMAX= atoi(argv[a]);
NW= NMAX+1;
LT= calloc(NW,sizeof(cnt_t));
GF= calloc(NW,sizeof(cnt_t));
EQ= calloc(NW,sizeof(cnt_t));

/*:85*/
#line 1911 "diag180.w"

clock_t t0= clock();
T_START= T_LAST= wall();
/*86:*/
#line 1970 "diag180.w"

for(int str= 0;str<2;str++){
FILE*f= str?stdout:stderr;
fprintf(f,"# command:");
for(int q= 0;q<argc;q++)
fprintf(f," %s",argv[q]);
fprintf(f,"\n# diag180 built %s %s; N=%d",__DATE__,__TIME__,NMAX);
if(argc-a==4)
fprintf(f,", single case c=%s d=%s vhi=%s",argv[a+1],argv[a+2],argv[a+3]);
if(MOD)
fprintf(f,", modulus %llu\n",(unsigned long long)MOD);
else
fprintf(f,", modulus 2**"CNT_BITS"\n");
fprintf(f,"# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
OLDBOUND?"old Kruskal (-K)":"cactus",HEIGHT?"on":"off (-H)",
LOOKAHEAD==2?"every cell (-A)":LOOKAHEAD?"at the kink":"off (-N)",
REFLECT?"on":"off (-R)",CHECK?"on (-C)":"off");
fprintf(f,"# output: results file %s, progress -v %d -p %g s\n",
RESFILE,VERB,PSEC);
fprintf(f,"# sizes: counter %zu bytes, count record %zu bytes, key word 64 bits\n",
sizeof(cnt_t),sizeof(rec_t));
fflush(f);
}

/*:86*/
#line 1914 "diag180.w"

if(argc-a==4){
/*87:*/
#line 1996 "diag180.w"

int c= atoi(argv[a+1]),d= atoi(argv[a+2]),vhi= atoi(argv[a+3]);
count_case(c,d,vhi);
if(MOD)
printf("# modulus %llu; n #(vspan<uspan) #(vspan=uspan)\n",(unsigned long long)MOD);
else
printf("# modulus 2**"CNT_BITS"; n #(vspan<uspan) #(vspan=uspan)\n");
for(int n= 1;n<=NMAX;n++)
printf("%d %s %s\n",n,cs(LT[n]),cs(EQ[n]));

/*:87*/
#line 1916 "diag180.w"

}else{
/*88:*/
#line 2009 "diag180.w"


cnt_t*out= calloc(3*NW,sizeof(cnt_t));
for(int c= 0;c<=1;c++)
for(int d= 0;d<=1;d++){



if(c==0&&d==1)
continue;
int t= c-d,ty= c+d,twice= c!=d;
int vhi= t>=0?(t+1)/2:0;
for(;2*vhi-t+1<=NMAX;vhi++){
if(2*vhi-t+1<1)
continue;
count_case(c,d,vhi);
for(int n= 1;n<=NMAX;n++){
cnt_t*o= &out[ty*NW+n];
*o= addm(*o,addm(LT[n],LT[n]));
*o= addm(*o,EQ[n]);
if(twice){
*o= addm(*o,addm(LT[n],LT[n]));
*o= addm(*o,EQ[n]);
}
}
}
}

/*:88*/
#line 1918 "diag180.w"

/*89:*/
#line 2045 "diag180.w"



FILE*rf= fopen(RESFILE,"a");
if(!rf){
perror(RESFILE);
return 1;
}
char mstr[32];
if(MOD)
snprintf(mstr,sizeof mstr,"%llu",(unsigned long long)MOD);
else
snprintf(mstr,sizeof mstr,"2**"CNT_BITS);
fprintf(rf,"# diag180 N=%d modulus %s\ntry:\n    Results\nexcept NameError:\n    Results = {}\n",NMAX,mstr);
static const char*tyname[3]= {"R180C","R180M","R180V"};
for(int n= 1;n<=NMAX;n++)
for(int ty= 0;ty<3;ty++)
fprintf(rf,"Results[%d, %s, '%s'] = %s\n",n,mstr,tyname[ty],
cs(out[ty*NW+n]));
fclose(rf);
printf("# modulus %s; results appended to %s\n",mstr,RESFILE);

/*:89*/
#line 1919 "diag180.w"

/*90:*/
#line 2069 "diag180.w"

printf("# n fix(r180) R180C(centre in a cell) R180M(edge midpoint) R180V(vertex)\n");
for(int n= 1;n<=NMAX;n++){
cnt_t tot= addm(addm(out[n],out[NW+n]),out[2*NW+n]);
printf("%d %s %s %s %s\n",n,cs(tot),cs(out[n]),cs(out[NW+n]),
cs(out[2*NW+n]));
}

/*:90*/
#line 1920 "diag180.w"

}
/*91:*/
#line 2077 "diag180.w"

fprintf(stderr,"# peak memory %.1f MB (both tables; %zu count records, key %zu 64-bit words)\n",
MAXBYTES/1e6,MAXRECS,MAXKW);
fprintf(stderr,"# max states %zu, %.2fs\n",MAXSTATES,
(double)(clock()-t0)/CLOCKS_PER_SEC);

/*:91*/
#line 1922 "diag180.w"

return 0;
}

/*:84*/
#line 109 "diag180.w"


/*:2*/
