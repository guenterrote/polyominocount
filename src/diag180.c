#define MAXW 160
#define PSELF 0xFF
#define PGHOST 0xFE
#define MAXLAB (MAXW+2) 
/*2:*/
#line 154 "diag180.w"

/*3:*/
#line 161 "diag180.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <stdint.h> 
#include <time.h> 

/*:3*/
#line 155 "diag180.w"

/*18:*/
#line 665 "diag180.w"

enum{OK,DEAD,COMPLETE};

/*:18*//*23:*/
#line 805 "diag180.w"

enum{C_EMPTY,C_CONTINUE,C_SINGLETON,C_FIRST,C_MIDDLE,C_LAST};

/*:23*//*28:*/
#line 929 "diag180.w"

#ifdef CNT128
typedef unsigned __int128 cnt_t;
#define CNT_BITS "128"
#else
#line 934 "diag180.w"
 typedef uint64_t cnt_t;
#define CNT_BITS "64"
#endif
#line 937 "diag180.w"

/*:28*//*32:*/
#line 990 "diag180.w"

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

/*:32*/
#line 156 "diag180.w"

/*8:*/
#line 340 "diag180.w"

static int NMAX;
static int M_x,M_y;
static int M_u,M_v;
static int umin,umax;
static int W;
static int E;
static int NMAX_plus_1;

/*:8*//*13:*/
#line 532 "diag180.w"

static int MAXC;
static int KL;

/*:13*//*22:*/
#line 802 "diag180.w"

static int KW;

/*:22*//*29:*/
#line 938 "diag180.w"

static uint64_t MOD= 0;

/*:29*//*48:*/
#line 1402 "diag180.w"

static int CUR_V,VSTAR,V0;
static int HEIGHT= 1;
static int H0,S0;
static int NOPT,OPTD[2*MAXW+4],OPTP[2*MAXW+4],OPTS[2*MAXW+4];
static int NCL,CLD[MAXW],CLP[MAXW],CLS[MAXW];
static int DEDGE;
static int CUR_TOUCHES;

/*:48*//*65:*/
#line 1889 "diag180.w"

static cnt_t*PORTRAIT,*SQUARE;
static cnt_t*RECT;
static size_t MAXSTATES,MAXBYTES,MAXRECS,MAXKW;

/*:65*//*69:*/
#line 1969 "diag180.w"

static int VERB= 1;
static double PSEC= 60;
static double T_START,T_LAST;

/*:69*//*81:*/
#line 2255 "diag180.w"

static int REFLECT= 1;
static int OLDBOUND= 0;
static const char*RESFILE= "results180.py";
static int LOOKAHEAD= 1;
static int CHECK= 0;

/*:81*/
#line 157 "diag180.w"

/*9:*/
#line 368 "diag180.w"

static int vspan(int k){return 2*k+1-E;}

/*:9*//*14:*/
#line 547 "diag180.w"

static inline int is_class(int p){return p!=PSELF&&p!=PGHOST;}

static inline void relabel(uint8_t*f,int from,int to)
{
for(int j= 0;j<W;j++)
if(f[j]==from)
f[j]= to;
}

/*:14*//*15:*/
#line 571 "diag180.w"

static void join(uint8_t*f,uint8_t*par,int a,int b)
{
if(a==b)
return;
int pa= par[a],pb= par[b];
relabel(f,b,a);
if(pa==b)
par[a]= PSELF;
else if(pa==PSELF||pb==PSELF){
if(is_class(pa))
relabel(f,pa,a);
if(is_class(pb))
relabel(f,pb,a);
par[a]= PSELF;
}else if(is_class(pa)){
if(is_class(pb))
relabel(f,pb,pa);
}else if(is_class(pb)){
par[a]= pb;
par[pb]= a;
}
}

/*:15*//*16:*/
#line 604 "diag180.w"

static void join_image(uint8_t*f,uint8_t*par,int a,int b)
{
int pa= par[a],pb= par[b];
if(pb==PSELF)
join(f,par,a,b);
else if(is_class(pb))
join(f,par,a,pb);
else if(a==b)
par[a]= PSELF;
else if(pa==PSELF)
relabel(f,b,a);
else if(is_class(pa))
relabel(f,b,pa);
else{
par[a]= b;
par[b]= a;
}
}

/*:16*//*17:*/
#line 626 "diag180.w"

static int canon(const uint8_t*f,const uint8_t*par,int touches,uint8_t*key)
{
uint8_t map[MAXLAB],old[MAXLAB];
int m= 0;
memset(map,0,sizeof(map));
memset(key,0,KL);
for(int i= 0;i<W;i++){
if(!f[i])
continue;
if(!map[f[i]]){
map[f[i]]= ++m;
old[m]= f[i];
}
key[i]= map[f[i]];
}
for(int j= 1;j<=m;j++){
int q= par[old[j]];
key[W+j-1]= is_class(q)?map[q]:q;
}
key[KL-1]= touches;
return m;
}

/*:17*//*19:*/
#line 668 "diag180.w"

static int retire(uint8_t*f,uint8_t*par,int i)
{
int lab= f[i];
if(!lab)
return OK;
f[i]= 0;
int any= 0;
for(int j= 0;j<W;j++)
if(f[j]){
any= 1;
if(f[j]==lab)
return OK;
}
int p= par[lab];
if(p==PSELF)
return any?DEAD:COMPLETE;
if(p==PGHOST)
return DEAD;
par[p]= PGHOST;
return OK;
}

/*:19*//*20:*/
#line 708 "diag180.w"

static int step(const uint8_t*sk,int touches,int i,int occ,int last,
int first_twin,uint8_t*key,int*result_touches,int*mp)
{
uint8_t f[MAXW],par[MAXLAB];
memcpy(f,sk,W);
memcpy(par+1,sk+W,MAXC);
if(occ){
int nw= 1;
while(nw<=MAXC&&par[nw])
nw++;
par[nw]= PGHOST;
f[i]= nw;
for(int nb= i-1;nb<=i+1;nb+= 2)
if(nb>=0&&nb<W&&f[nb])
join(f,par,nw,f[nb]);
if(first_twin){
for(int jj= W-2-i;jj<=W-i;jj+= 2)
if(jj>=0&&jj<=i&&f[jj])
join_image(f,par,nw,f[jj]);
}
}
int status= OK;
if(i-1>=0)
status= retire(f,par,i-1);
if(status==OK&&last&&i+1<W)
status= retire(f,par,i+1);
if(status!=OK)
return status;
*result_touches= touches||(occ&&(i==0||i==W-1));
*mp= canon(f,par,*result_touches,key);
return OK;
}

/*:20*//*25:*/
#line 829 "diag180.w"

static inline void putbits(uint64_t*w,int*pos,uint64_t val,int nb)
{
int i= *pos>>6,o= *pos&63;
w[i]|= val<<o;
if(o+nb> 64)
w[i+1]|= val>>(64-o);
*pos+= nb;
}

static inline int getbits(const uint64_t*w,int*pos,int nb)
{
int i= *pos>>6,o= *pos&63;
uint64_t val= w[i]>>o;
if(o+nb> 64)
val|= w[i+1]<<(64-o);
*pos+= nb;
return val&((1u<<nb)-1);
}

/*:25*//*26:*/
#line 851 "diag180.w"

static void pack(const uint8_t*key,uint64_t*w)
{
int first[MAXW+2],lastrun[MAXW+2];
int pos= 0,m= 0;
memset(w,0,KW*8);
for(int i= 0;i<W;i++)
if(key[i]&&(i==0||!key[i-1])){
int l= key[i];
if(l> m){
m= l;
first[l]= i;
}
lastrun[l]= i;
}
for(int i= 0;i<W;i++){
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
putbits(w,&pos,key[W+j-1]!=PGHOST,1);
putbits(w,&pos,key[KL-1],1);
}

/*:26*//*27:*/
#line 885 "diag180.w"

static void unpack(const uint64_t*w,uint8_t*key)
{
int stack[MAXW+2],sp= 0,m= 0,pos= 0;
memset(key,0,KL);
for(int i= 0;i<W;i++){
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
key[W+j-1]= PGHOST;
for(int a= 0,b= sp-1;a<=b;a++,b--){
key[W+stack[a]-1]= a==b?PSELF:stack[b];
key[W+stack[b]-1]= a==b?PSELF:stack[a];
}
key[KL-1]= getbits(w,&pos,1);
}

/*:27*//*30:*/
#line 941 "diag180.w"

static inline cnt_t addm(cnt_t a,cnt_t b)
{
if(!MOD)
return a+b;
cnt_t s= a+b;
if(s<a||s>=MOD)
s-= MOD;
return s;
}

/*:30*//*31:*/
#line 956 "diag180.w"

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

/*:31*//*33:*/
#line 1006 "diag180.w"

static void*xrealloc(void*p,size_t s)
{
p= realloc(p,s);
if(!p){
fprintf(stderr,"out of memory\n");
exit(1);
}
return p;
}

/*:33*//*34:*/
#line 1017 "diag180.w"

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

/*:34*//*35:*/
#line 1057 "diag180.w"

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

/*:35*//*36:*/
#line 1071 "diag180.w"

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

/*:36*//*37:*/
#line 1104 "diag180.w"

static size_t tab_get(table_t*t,const uint8_t*key)
{
uint64_t w[12];
pack(key,w);
if(CHECK){
uint8_t k2[2*MAXW+4];
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

/*:37*//*38:*/
#line 1142 "diag180.w"

static void tab_add(table_t*t,size_t si,int wt,cnt_t c)
{
uint32_t prev= 0,r= t->head[si];
while(r&&t->rec[r-1].wt<wt){
prev= r;
r= t->rec[r-1].next;
}
if(r&&t->rec[r-1].wt==wt){
t->rec[r-1].c= addm(t->rec[r-1].c,c);
return;
}
if(t->nrec==t->rcap){
if(t->rcap>=0xFFFFFFFFu){
fprintf(stderr,"too many count records for 32-bit links\n");
exit(1);
}
t->rcap= 2*t->rcap> 0xFFFFFFFFu?0xFFFFFFFFu:2*t->rcap;
t->rec= xrealloc(t->rec,t->rcap*sizeof(rec_t));
}
rec_t*q= &t->rec[t->nrec++];
q->c= c;
q->wt= wt;
q->next= r;
if(prev)
t->rec[prev-1].next= t->nrec;
else
t->head[si]= t->nrec;
}

/*:38*//*39:*/
#line 1180 "diag180.w"

static void insert_pruned(table_t*t,const uint8_t*key,const table_t*src,
size_t si,int add,int lowerbound)
{
int lim= NMAX-lowerbound-add;
uint32_t r= src->head[si];
if(!r||src->rec[r-1].wt> lim)
return;
size_t ti= tab_get(t,key);
for(;r&&src->rec[r-1].wt<=lim;r= src->rec[r-1].next)
tab_add(t,ti,src->rec[r-1].wt+add,src->rec[r-1].c);
}

/*:39*//*41:*/
#line 1208 "diag180.w"

static int completion_possible(const uint8_t*key)
{
int m= 0;
for(int j= 0;j<MAXC;j++){
int q= key[W+j];
if(!q)
break;
m++;
if(q!=PGHOST)
return 1;
}
return m==0;
}

/*:41*//*49:*/
#line 1413 "diag180.w"

static inline int frow(int j)
{
return((umin+j-CUR_V)%2+2)%2?CUR_V-1:CUR_V;
}

/*:49*//*50:*/
#line 1419 "diag180.w"

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

/*:50*//*51:*/
#line 1434 "diag180.w"

static inline int both_cost(int h,int j)
{
int e= j<W-1-j?j:W-1-j;
int r= VSTAR-h;
return e> r?e:r;
}

/*:51*//*52:*/
#line 1445 "diag180.w"

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

/*:52*//*53:*/
#line 1464 "diag180.w"

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

/*:53*//*54:*/
#line 1507 "diag180.w"

static inline int finish(int h,int s)
{
int top= VSTAR-h> 0?VSTAR-h:0;
if(CUR_TOUCHES)
return HEIGHT?top:0;
if(!HEIGHT)
return DEDGE;
int sep= DEDGE+top;
return s<sep?s:sep;
}

/*:54*//*55:*/
#line 1521 "diag180.w"

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

/*:55*//*57:*/
#line 1620 "diag180.w"

static int cactus_bound(const uint8_t*key,int touches,int kink)
{
/*58:*/
#line 1634 "diag180.w"

int rs[MAXW],re[MAXW],rl[MAXW],R= 0;
int lastrun[MAXW+2];
for(int i= 0;i<W;i++)
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

/*:58*/
#line 1623 "diag180.w"

/*59:*/
#line 1650 "diag180.w"

H0= -1000;
S0= 1<<20;
NOPT= 0;
NCL= 0;
CUR_TOUCHES= touches;
for(int i= 0;i<W;i++)
if(key[i]){
reach(frow(i));
sreach(both_cost(frow(i),i));
}

/*:59*/
#line 1624 "diag180.w"

int cost= 0;
/*60:*/
#line 1669 "diag180.w"


int sl[MAXW],sb[MAXW],sp= 0;
int gb[MAXW],pb[MAXW],qb[MAXW],nb= 0;
int toplab[MAXW],topgap[MAXW],toppk[MAXW],topsc[MAXW],nt= 0;
for(int j= 0;j<R;j++){
int X= rl[j];
if(j> 0){
int a= re[j-1],b= rs[j];
int g= b-a-1+(b==kink);

int pk= (g-1+(b==kink)+frow(a)+frow(b)+2)/2;
int sc= touches?1<<20:gap_score(a,frow(a),b==kink?b+1:b-1,frow(b)+1);
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

/*:60*/
#line 1626 "diag180.w"

/*61:*/
#line 1708 "diag180.w"


int cr[MAXW],m= 0;
for(int t= 0;t<nt;t++)
if(key[W+toplab[t]-1]!=PGHOST)
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
int g[2*MAXW],pk[2*MAXW],scs[2*MAXW],n= 0,m1,p1,s1;
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

/*:61*/
#line 1627 "diag180.w"

/*62:*/
#line 1761 "diag180.w"

DEDGE= 0;
if(!touches){
int d1= rs[0],d2= W-1-re[R-1];
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

/*:62*/
#line 1628 "diag180.w"

return cost;
}

/*:57*//*63:*/
#line 1784 "diag180.w"

static int exact_bound(const uint8_t*key,int touches)
{
int sup[MAXW+1];
int occ[MAXW],no= 0;
for(int j= 1;j<=MAXC;j++){
int q= key[W+j-1];
if(!q)
break;
sup[j]= (q!=PSELF&&q!=PGHOST&&q<j)?q:j;
}
for(int i= 0;i<W;i++)
if(key[i])
occ[no++]= i;
if(!no)
return 0;

int par[MAXW+1];
int ncls= 0;
int seen[MAXW+1];
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
int wgt[MAXW],ord[MAXW],bucket[MAXW+1];
memset(bucket,0,sizeof(bucket));
for(int a= 0;a<ne;a++){
wgt[a]= occ[a+1]-occ[a]-1;
bucket[wgt[a]]++;
}
for(int g= 1;g<=W;g++)
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
if(!touches){
int d1= occ[0],d2= W-1-occ[no-1];
cost+= d1<d2?d1:d2;
}
return cost;
}

/*:63*//*64:*/
#line 1855 "diag180.w"

static int bound(const uint8_t*key,int touches,int k,int kmin,int i,int last,int first_twin)
{
int lowerbound= kmin> k?2*(kmin-k):0;
if(!touches){
int dist= W;
for(int j= 0;j<W;j++)
if(key[j]){
int dd= j<W-1-j?j:W-1-j;
if(dd<dist)
dist= dd;
}
if(2*dist> lowerbound)
lowerbound= 2*dist;
}
if(!first_twin||last){
if(!completion_possible(key))
return-1;
int kink= (!last&&!key[i]&&i+1<W&&key[i+1])?i+1:-1;
CUR_V= V0+k-1;
VSTAR= V0+kmin-1;
int eb= 2*(OLDBOUND?exact_bound(key,touches):cactus_bound(key,touches,kink));
if(eb> lowerbound)
lowerbound= eb;
}
return lowerbound;
}

/*:64*//*67:*/
#line 1903 "diag180.w"

static void record(int ktop,const table_t*src,size_t si,int touches)
{
int h= vspan(ktop);
if(!touches||h<W)
return;
cnt_t*tgt= h==W?SQUARE:PORTRAIT;
for(uint32_t r= src->head[si];r;r= src->rec[r-1].next){
int wt= src->rec[r-1].wt;
tgt[wt]= addm(tgt[wt],src->rec[r-1].c);
RECT[wt]= addm(RECT[wt],src->rec[r-1].c);
}
}

/*:67*//*68:*/
#line 1942 "diag180.w"

static void rect_line(int ktop)
{
static const char*tyname[3]= {
"R180C center in a cell",
"R180M center on the midpoint of a vertical edge",
"R180V center at a vertex"};
int h= vspan(ktop);
int vtop= E?(M_v+1)/2+ktop-1:M_v/2+ktop;
if(VERB>=1){
fprintf(stderr,"#rect %s (Mx=%d My=%d): width %d (u=%d..%d) x height %d (v=%d..%d):",
tyname[M_x+M_y],M_x,M_y,W,umin,umin+W-1,h,M_v-vtop,vtop);
for(int n= 0;n<=NMAX;n++)
if(RECT[n])
fprintf(stderr," n=%d:%s",n,cs(RECT[n]));
fprintf(stderr,"\n");
fflush(stderr);
}
memset(RECT,0,NMAX_plus_1*sizeof(cnt_t));
}

/*:68*//*70:*/
#line 1974 "diag180.w"

static double wall(void)
{
struct timespec ts;
clock_gettime(CLOCK_MONOTONIC,&ts);
return ts.tv_sec+1e-9*ts.tv_nsec;
}
static void stat_line(const char*tag,int k,int v,int ci,int nu,
const table_t*cur,const table_t*nxt)
{
double now= wall();
T_LAST= now;
cnt_t done= 0;
for(int n= 0;n<=NMAX;n++)
done= addm(done,addm(PORTRAIT[n],SQUARE[n]));
fprintf(stderr,"%s Mx=%d My=%d umax=%d W=%d row k=%d v=%d cell %d/%d: "
"%zu states, %zu records, %.1f MB, done %s, %.1fs wall, %.1fs cpu\n",
tag,M_x,M_y,umax,W,k,v,ci,nu,cur->n,cur->nrec,
(tab_bytes(cur)+tab_bytes(nxt))/1e6,cs(done),
now-T_START,(double)clock()/CLOCKS_PER_SEC);
fflush(stderr);
}

/*:70*//*71:*/
#line 1998 "diag180.w"

static void count_case(void)
{
/*72:*/
#line 2021 "diag180.w"

M_v= M_x+M_y;
M_u= M_x-M_y;
umin= M_u-umax;
W= umax-umin+1;
E= M_v&1;
memset(PORTRAIT,0,NMAX_plus_1*sizeof(cnt_t));
memset(SQUARE,0,NMAX_plus_1*sizeof(cnt_t));
if(W==1){
if(M_x==0&&M_y==0){
SQUARE[1]= addm(SQUARE[1],1);
RECT[1]= 1;
rect_line(0);
}
return;
}
if(W> MAXW){
fprintf(stderr,"strip too wide\n");
exit(1);
}
MAXC= (W+1)/2+1;
KL= W+MAXC+1;
KW= (3*W+MAXC+1+63)/64;
NMAX_plus_1= NMAX+1;

int kmin= 0;
while(vspan(kmin)<W)
kmin++;

/*:72*/
#line 2001 "diag180.w"

/*73:*/
#line 2052 "diag180.w"

table_t A,B,*cur= &A,*nxt= &B;
tab_init(&A,1024);
tab_init(&B,1024);
uint8_t key[MAXW*2+4];
int v0;

/*:73*/
#line 2002 "diag180.w"

/*74:*/
#line 2066 "diag180.w"

if(E==0){
int vm= M_v/2;
int half[MAXW],nh= 0;
for(int i= 0;i<W;i++)
if(((umin+i-vm)%2+2)%2==0&&2*(umin+i)>=M_u)
half[nh++]= i;

for(uint64_t mask= 1;mask<(1ull<<nh);mask++){
int cells[MAXW],nc= 0,isc[MAXW];
memset(isc,0,sizeof(int)*W);
for(int b= 0;b<nh;b++)
if(mask>>b&1){
isc[half[b]]= 1;
isc[W-1-half[b]]= 1;
}
for(int i= 0;i<W;i++)
if(isc[i])
cells[nc++]= i;
if(nc> NMAX)
continue;
uint8_t f[MAXW],par[MAXLAB];
memset(f,0,W);
for(int a= 0;a<nc;a++)
f[cells[a]]= a+1;
for(int a= 0;a<nc;a++){
int i= cells[a],j= W-1-i;
par[a+1]= j==i?PSELF:f[j];
}
int touches= isc[0]||isc[W-1];
canon(f,par,touches,key);
tab_add(cur,tab_get(cur,key),nc,1);
}
v0= vm+1;
V0= v0;
}else{
memset(key,0,KL);
tab_add(cur,tab_get(cur,key),0,1);
v0= (M_v+1)/2;
V0= v0;
}

/*:74*/
#line 2003 "diag180.w"

int cells_v[MAXW];
for(int k= 1;cur->n;k++){
/*75:*/
#line 2110 "diag180.w"

int v= v0+k-1,nu= 0;
for(int i= 0;i<W;i++)
if(((umin+i-v)%2+2)%2==0)
cells_v[nu++]= i;
int first_twin= (E==1&&k==1);

/*:75*/
#line 2006 "diag180.w"

for(int ci= 0;ci<nu;ci++){
/*76:*/
#line 2117 "diag180.w"

int i= cells_v[ci];
int last= ci==nu-1;
tab_clear(nxt);
for(size_t si= 0;si<cur->n;si++){
/*77:*/
#line 2128 "diag180.w"

uint8_t sk[MAXW*2+4];
unpack(cur->keys+si*KW,sk);
int touches= sk[KL-1];
int empty= 1;
for(int j= 0;j<W;j++)
if(sk[j]){
empty= 0;
break;
}
for(int occ= 0;occ<=1;occ++){
if(occ&&k> 1&&empty)
continue;
int result_touches,m;
int status= step(sk,touches,i,occ,last,first_twin,key,&result_touches,&m);
if(status==DEAD)
continue;
if(status==COMPLETE){
record(k-1,cur,si,touches);
continue;
}
int add= occ?2:0;
if(m==0){
if(first_twin&&!last)
insert_pruned(nxt,key,cur,si,add,0);
continue;
}
int lowerbound= bound(key,result_touches,k,kmin,i,last,first_twin);
if(lowerbound<0)
continue;
/*78:*/
#line 2172 "diag180.w"

if(LOOKAHEAD&&!last&&!first_twin&&(LOOKAHEAD==2||(!key[i]&&i+1<W&&key[i+1]))){

int i2= cells_v[ci+1],last2= ci+1==nu-1;
int la= 1<<20;
uint8_t k2[MAXW*2+4];
for(int o2= 0;o2<=1&&la;o2++){
int result_touches2,m2;
int st2= step(key,result_touches,i2,o2,last2,0,k2,&result_touches2,&m2);
if(st2==DEAD)
continue;

if(st2==COMPLETE&&!(result_touches&&vspan(k-1)>=W))
continue;
int val= 2*o2;
if(st2==OK){
if(m2==0)
continue;
int lowerbound2= bound(k2,result_touches2,k,kmin,i2,last2,0);
if(lowerbound2<0)
continue;
val+= lowerbound2;
}
if(val<la)
la= val;
}
if(la> lowerbound)
lowerbound= la;
}

/*:78*/
#line 2158 "diag180.w"

insert_pruned(nxt,key,cur,si,add,lowerbound);
}

/*:77*/
#line 2122 "diag180.w"

}
/*79:*/
#line 2202 "diag180.w"

table_t*x= cur;
cur= nxt;
nxt= x;
if(cur->n> MAXSTATES)
MAXSTATES= cur->n;
size_t by= tab_bytes(cur)+tab_bytes(nxt);
if(by> MAXBYTES)
MAXBYTES= by,MAXRECS= cur->nrec,MAXKW= KW;
if(!last&&(VERB>=2||(VERB==1&&wall()-T_LAST>=PSEC)))
stat_line("#cell",k,v,ci+1,nu,cur,nxt);

/*:79*/
#line 2124 "diag180.w"


/*:76*/
#line 2008 "diag180.w"

}
/*80:*/
#line 2222 "diag180.w"

if(REFLECT&&M_u==0&&cur->n){

tab_clear(nxt);
uint8_t rkey[MAXW*2+4];
for(size_t si= 0;si<cur->n;si++){
uint8_t sk[MAXW*2+4];
unpack(cur->keys+si*KW,sk);
uint8_t f[MAXW];
for(int a= 0;a<W;a++)
f[a]= sk[W-1-a];
canon(f,sk+W-1,sk[KL-1],rkey);
const uint8_t*use= memcmp(rkey,sk,KL)<0?rkey:sk;
insert_pruned(nxt,use,cur,si,0,0);
}
table_t*x= cur;
cur= nxt;
nxt= x;
}

/*:80*/
#line 2010 "diag180.w"

if(VERB>=1)
stat_line("#row ",k,v,nu,nu,cur,nxt);
if(vspan(k-1)>=W)
rect_line(k-1);
}
tab_free(&A);
tab_free(&B);
}

/*:71*/
#line 158 "diag180.w"

/*82:*/
#line 2262 "diag180.w"

int main(int argc,char**argv)
{
/*83:*/
#line 2280 "diag180.w"

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
fprintf(stderr,"usage: %s [options] N [Mx My umax]\n",argv[0]);
return 1;
}
NMAX= atoi(argv[a]);
if(NMAX<1||NMAX> MAXW){
fprintf(stderr,"%s: N must be between 1 and %d\n",argv[0],MAXW);
return 1;
}
NMAX_plus_1= NMAX+1;
/*66:*/
#line 1896 "diag180.w"

PORTRAIT= calloc(NMAX_plus_1,sizeof(cnt_t));
SQUARE= calloc(NMAX_plus_1,sizeof(cnt_t));
RECT= calloc(NMAX_plus_1,sizeof(cnt_t));

/*:66*/
#line 2320 "diag180.w"


/*:83*/
#line 2265 "diag180.w"

clock_t t0= clock();
T_START= T_LAST= wall();
/*85:*/
#line 2350 "diag180.w"

for(int str= 0;str<2;str++){
FILE*f= str?stdout:stderr;
fprintf(f,"# command:");
for(int q= 0;q<argc;q++)
fprintf(f," %s",argv[q]);
fprintf(f,"\n# diag180 built %s %s; N=%d",__DATE__,__TIME__,NMAX);
if(argc-a==4)
fprintf(f,", single case Mx=%s My=%s umax=%s",argv[a+1],argv[a+2],argv[a+3]);
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

/*:85*/
#line 2268 "diag180.w"

if(argc-a==4){
/*86:*/
#line 2376 "diag180.w"

M_x= atoi(argv[a+1]),M_y= atoi(argv[a+2]),umax= atoi(argv[a+3]);
count_case();
if(MOD)
printf("# modulus %llu; n #(uspan<vspan) #(uspan=vspan)\n",(unsigned long long)MOD);
else
printf("# modulus 2**"CNT_BITS"; n #(uspan<vspan) #(uspan=vspan)\n");
for(int n= 1;n<=NMAX;n++)
printf("%d %s %s\n",n,cs(PORTRAIT[n]),cs(SQUARE[n]));

/*:86*/
#line 2270 "diag180.w"

}else{
/*87:*/
#line 2389 "diag180.w"


cnt_t*out= calloc(3*NMAX_plus_1,sizeof(cnt_t));
for(M_x= 0;M_x<=1;M_x++)
for(M_y= 0;M_y<=M_x;M_y++){


int t= M_x-M_y,ty= M_x+M_y,twice= M_x!=M_y;
for(umax= t;2*umax-t+1<=NMAX;umax++){
count_case();
for(int n= 1;n<=NMAX;n++){
cnt_t*o= &out[ty*NMAX_plus_1+n];
*o= addm(*o,addm(PORTRAIT[n],PORTRAIT[n]));
*o= addm(*o,SQUARE[n]);
if(twice){
*o= addm(*o,addm(PORTRAIT[n],PORTRAIT[n]));
*o= addm(*o,SQUARE[n]);
}
}
}
}

/*:87*/
#line 2272 "diag180.w"

/*88:*/
#line 2419 "diag180.w"



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
cs(out[ty*NMAX_plus_1+n]));
fclose(rf);
printf("# modulus %s; results appended to %s\n",mstr,RESFILE);

/*:88*/
#line 2273 "diag180.w"

/*89:*/
#line 2443 "diag180.w"

printf("# n fix(r180) R180C(center in a cell) R180M(edge midpoint) R180V(vertex)\n");
for(int n= 1;n<=NMAX;n++){
cnt_t tot= addm(addm(out[n],out[NMAX_plus_1+n]),out[2*NMAX_plus_1+n]);
printf("%d %s %s %s %s\n",n,cs(tot),cs(out[n]),cs(out[NMAX_plus_1+n]),
cs(out[2*NMAX_plus_1+n]));
}

/*:89*/
#line 2274 "diag180.w"

}
/*90:*/
#line 2451 "diag180.w"

fprintf(stderr,"# peak memory %.1f MB (both tables; %zu count records, key %zu 64-bit words)\n",
MAXBYTES/1e6,MAXRECS,MAXKW);
fprintf(stderr,"# max states %zu, %.2fs\n",MAXSTATES,
(double)(clock()-t0)/CLOCKS_PER_SEC);

/*:90*/
#line 2276 "diag180.w"

return 0;
}

/*:82*/
#line 159 "diag180.w"


/*:2*/
