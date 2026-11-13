#define max(a,b) ((a) > (b) ?(a) :(b) ) 
#define min(a,b) ((a) <(b) ?(a) :(b) )  \
 \

#define MAXW 160
#define MAXLAB (MAXW+2) 
#define NO_STATE ((size_t) -1) 
/*2:*/
#line 167 "diag180.w"

/*3:*/
#line 174 "diag180.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <stdint.h> 
#include <time.h> 

/*:3*/
#line 168 "diag180.w"

/*13:*/
#line 646 "diag180.w"

typedef struct{
uint8_t label[MAXW];
uint8_t partner[MAXLAB];
uint8_t num_classes;
uint8_t touches;
}state_t;

/*:13*//*20:*/
#line 794 "diag180.w"

enum{OK,DEAD,COMPLETE};

/*:20*//*25:*/
#line 928 "diag180.w"

enum{C_EMPTY,C_CONTINUE,C_SINGLETON,C_FIRST,C_MIDDLE,C_LAST};

/*:25*//*28:*/
#line 978 "diag180.w"

typedef struct{
uint8_t code[MAXW];
uint8_t crossing[MAXLAB];
uint8_t num_classes;
uint8_t touches;
}coded_t;

/*:28*//*33:*/
#line 1109 "diag180.w"

#ifdef CNT128
typedef unsigned __int128 cnt_t;
#define CNT_BITS "128"
#else
#line 1114 "diag180.w"
 typedef uint64_t cnt_t;
#define CNT_BITS "64"
#endif
#line 1117 "diag180.w"

/*:33*//*37:*/
#line 1174 "diag180.w"

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
uint8_t*lowerbound;
}table_t;

/*:37*/
#line 169 "diag180.w"

/*8:*/
#line 362 "diag180.w"

static int NMAX;
static int M_x,M_y;
static int M_u,M_v;
static int umin,umax;
static int W;
static int E;
static int NMAX_plus_1;

/*:8*//*14:*/
#line 654 "diag180.w"

static int MAXC;

/*:14*//*24:*/
#line 925 "diag180.w"

static int KW;

/*:24*//*34:*/
#line 1118 "diag180.w"

static uint64_t MOD= 0;

/*:34*//*56:*/
#line 1645 "diag180.w"

static int CUR_V,VSTAR,V0;
#line 1648 "diag180.w"
static int H0,S0;
static int NOPT,OPTD[2*MAXW+4],OPTP[2*MAXW+4],OPTS[2*MAXW+4];
static int NCL,CLD[MAXW],CLP[MAXW],CLS[MAXW];
static int DEDGE;
static int CUR_TOUCHES;

/*:56*//*72:*/
#line 2256 "diag180.w"

static cnt_t*PORTRAIT,*SQUARE;
static cnt_t*RECT;
static size_t MAXSTATES,MAXBYTES,MAXRECS,MAXKW;

/*:72*//*76:*/
#line 2336 "diag180.w"

static int VERB= 1;
static double PSEC= 60;
static double T_START,T_LAST;

/*:76*//*89:*/
#line 2641 "diag180.w"

static int REFLECT= 1;
#line 2644 "diag180.w"
static const char*RESFILE= "results180.py";
#line 2646 "diag180.w"
static int CHECK= 0;

/*:89*/
#line 170 "diag180.w"

/*9:*/
#line 390 "diag180.w"

static int vspan(int k){return 2*k+1-E;}

/*:9*//*15:*/
#line 659 "diag180.w"

static int state_cmp(const state_t*a,const state_t*b)
{
int c= memcmp(a->label,b->label,W);
if(c)
return c;
c= memcmp(a->partner+1,b->partner+1,a->num_classes);
if(c)
return c;
return a->touches-b->touches;
}

/*:15*//*16:*/
#line 682 "diag180.w"

static inline int is_central(const uint8_t*par,int a){return par[a]==a;}
static inline int is_floating(const uint8_t*par,int a){return par[a]==0;}
static inline int is_paired(const uint8_t*par,int a){return par[a]!=0&&par[a]!=a;}

static inline void relabel(uint8_t*f,int from,int to)
{
for(int j= 0;j<W;j++)
if(f[j]==from)
f[j]= to;
}

/*:16*//*17:*/
#line 708 "diag180.w"

static void join(uint8_t*f,uint8_t*par,int a,int b)
{
if(a==b)
return;
int pa= par[a],pb= par[b];
relabel(f,b,a);
if(pa==b)
par[a]= a;
else if(pa==a||pb==b){
if(pa!=0&&pa!=a)
relabel(f,pa,a);
if(pb!=0&&pb!=b)
relabel(f,pb,a);
par[a]= a;
}else if(pa!=0){
if(pb!=0)
relabel(f,pb,pa);
}else if(pb!=0){
par[a]= pb;
par[pb]= a;
}
}

/*:17*//*18:*/
#line 741 "diag180.w"

static void join_image(uint8_t*f,uint8_t*par,int a,int b)
{
if(!is_floating(par,b))
join(f,par,a,par[b]);
else if(a==b)
par[a]= a;
else if(is_central(par,a))
relabel(f,b,a);
else if(is_paired(par,a))
relabel(f,b,par[a]);
else{
par[a]= b;
par[b]= a;
}
}

/*:18*//*19:*/
#line 760 "diag180.w"

static void canon(const uint8_t*f,const uint8_t*par,int touches,state_t*s)
{
uint8_t map[MAXLAB],old[MAXLAB];
int m= 0;
memset(map,0,sizeof(map));
for(int i= 0;i<W;i++){
if(f[i]&&!map[f[i]]){
map[f[i]]= ++m;
old[m]= f[i];
}
s->label[i]= map[f[i]];
}
for(int j= 1;j<=m;j++)
s->partner[j]= map[par[old[j]]];
s->num_classes= m;
s->touches= touches;
}

/*:19*//*21:*/
#line 797 "diag180.w"

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
if(is_central(par,lab))
return any?DEAD:COMPLETE;
if(is_floating(par,lab))
return DEAD;
par[par[lab]]= 0;
return OK;
}

/*:21*//*22:*/
#line 834 "diag180.w"

static int step(const state_t*s,int i,int occ,int last,int first_twin,
state_t*result)
{
uint8_t f[MAXW],par[MAXLAB];
memcpy(f,s->label,W);
memcpy(par,s->partner,s->num_classes+1);
if(occ){
int nw= s->num_classes+1;
par[nw]= 0;
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
canon(f,par,s->touches||(occ&&(i==0||i==W-1)),result);
return OK;
}

/*:22*//*27:*/
#line 952 "diag180.w"

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

/*:27*//*29:*/
#line 988 "diag180.w"

static void encode(const state_t*s,coded_t*c)
{
int first[MAXLAB],lastrun[MAXLAB];
int m= 0;
for(int i= 0;i<W;i++)
if(s->label[i]&&(i==0||!s->label[i-1])){
int l= s->label[i];
if(l> m){
m= l;
first[l]= i;
}
lastrun[l]= i;
}
for(int i= 0;i<W;i++){
int l= s->label[i];
if(!l)
c->code[i]= C_EMPTY;
else if(i> 0&&s->label[i-1])
c->code[i]= C_CONTINUE;
else if(first[l]==i)
c->code[i]= lastrun[l]==i?C_SINGLETON:C_FIRST;
else
c->code[i]= lastrun[l]==i?C_LAST:C_MIDDLE;
}
for(int j= 1;j<=m;j++)
c->crossing[j]= s->partner[j]!=0;
c->num_classes= m;
c->touches= s->touches;
}

/*:29*//*30:*/
#line 1023 "diag180.w"

static void decode(const coded_t*c,state_t*s)
{
int stack[MAXLAB],sp= 0,m= 0;
for(int i= 0;i<W;i++){
switch(c->code[i]){
case C_EMPTY:s->label[i]= 0;break;
case C_CONTINUE:s->label[i]= s->label[i-1];break;
case C_SINGLETON:s->label[i]= ++m;break;
case C_FIRST:s->label[i]= ++m;stack[sp++]= m;break;
case C_MIDDLE:s->label[i]= stack[sp-1];break;
case C_LAST:s->label[i]= stack[--sp];break;
}
}
sp= 0;
for(int j= 1;j<=m;j++)
if(c->crossing[j])
stack[sp++]= j;
else
s->partner[j]= 0;
for(int a= 0,b= sp-1;a<=b;a++,b--){
s->partner[stack[a]]= stack[b];
s->partner[stack[b]]= stack[a];
}
s->num_classes= m;
s->touches= c->touches;
}

/*:30*//*31:*/
#line 1053 "diag180.w"

static void compress(const coded_t*c,uint64_t*w)
{
int pos= 0;
memset(w,0,KW*8);
for(int i= 0;i<W;i++)
putbits(w,&pos,c->code[i],3);
for(int j= 1;j<=c->num_classes;j++)
putbits(w,&pos,c->crossing[j],1);
putbits(w,&pos,c->touches,1);
}

static void expand(const uint64_t*w,coded_t*c)
{
int pos= 0,m= 0;
for(int i= 0;i<W;i++){
c->code[i]= getbits(w,&pos,3);
if(c->code[i]==C_SINGLETON||c->code[i]==C_FIRST)
m++;
}
for(int j= 1;j<=m;j++)
c->crossing[j]= getbits(w,&pos,1);
c->num_classes= m;
c->touches= getbits(w,&pos,1);
}

/*:31*//*32:*/
#line 1079 "diag180.w"

static void pack(const state_t*s,uint64_t*w)
{
coded_t c;
encode(s,&c);
compress(&c,w);
}

static void unpack(const uint64_t*w,state_t*s)
{
coded_t c;
expand(w,&c);
decode(&c,s);
}

/*:32*//*35:*/
#line 1121 "diag180.w"

static inline cnt_t addm(cnt_t a,cnt_t b)
{
if(!MOD)
return a+b;
cnt_t s= a+b;
if(s<a||s>=MOD)
s-= MOD;
return s;
}

/*:35*//*36:*/
#line 1136 "diag180.w"

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

/*:36*//*38:*/
#line 1191 "diag180.w"

static void*xrealloc(void*p,size_t s)
{
p= realloc(p,s);
if(!p){
fprintf(stderr,"out of memory\n");
exit(1);
}
return p;
}

/*:38*//*39:*/
#line 1202 "diag180.w"

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
t->hash= 0;
t->lowerbound= 0;
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
free(t->lowerbound);
}

static size_t tab_bytes(const table_t*t)
{
return t->n*(KW*8+4)+t->nrec*sizeof(rec_t)
+(t->hash?(t->hmask+1)*4+t->n:0);
}

/*:39*//*40:*/
#line 1245 "diag180.w"

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

/*:40*//*41:*/
#line 1262 "diag180.w"

static void tab_grow(table_t*t)
{
size_t cap= t->cap*2;
t->keys= xrealloc(t->keys,cap*KW*8);
t->head= xrealloc(t->head,cap*4);
t->lowerbound= xrealloc(t->lowerbound,cap);
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


/*:41*//*42:*/
#line 1297 "diag180.w"

static void tab_clear(table_t*t)
{
t->n= 0;
t->nrec= 0;
t->hash= calloc(t->hmask+1,sizeof(uint32_t));
if(!t->hash){
fprintf(stderr,"out of memory\n");
exit(1);
}
t->lowerbound= xrealloc(0,t->cap);
}

static void tab_finish(table_t*t)
{
free(t->hash);
t->hash= 0;
free(t->lowerbound);
t->lowerbound= 0;
}

/*:42*//*43:*/
#line 1325 "diag180.w"

static size_t HASH_SLOT;

static size_t tab_find(const table_t*t,const uint64_t*w)
{
size_t p= hkey(w)&t->hmask;
while(t->hash[p]){
size_t i= t->hash[p]-1;
if(!memcmp(t->keys+i*KW,w,KW*8))
return i;
p= (p+1)&t->hmask;
}
HASH_SLOT= p;
return NO_STATE;
}

static size_t tab_create(table_t*t,const uint64_t*w)
{
size_t p= HASH_SLOT;
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
t->lowerbound[i]= 0;
return i;
}

/*:43*//*44:*/
#line 1360 "diag180.w"

static void pack_key(const state_t*s,uint64_t*w)
{
pack(s,w);
if(CHECK){
state_t s2;
unpack(w,&s2);
if(s2.num_classes!=s->num_classes||state_cmp(&s2,s)){
fprintf(stderr,"pack/unpack mismatch\n");
exit(2);
}
}
}

/*:44*//*45:*/
#line 1376 "diag180.w"

static size_t tab_get(table_t*t,const state_t*s)
{
uint64_t w[12];
pack_key(s,w);
size_t i= tab_find(t,w);
return i!=NO_STATE?i:tab_create(t,w);
}

/*:45*//*46:*/
#line 1388 "diag180.w"

static uint32_t new_record(table_t*t,int wt,cnt_t c,uint32_t next)
{
if(t->nrec==t->rcap){
if(t->rcap>=0xFFFFFFFFu){
fprintf(stderr,"too many count records for 32-bit links\n");
exit(1);
}
t->rcap= min(2*t->rcap,0xFFFFFFFFu);
t->rec= xrealloc(t->rec,t->rcap*sizeof(rec_t));
}
rec_t*q= &t->rec[t->nrec++];
q->c= c;
q->wt= wt;
q->next= next;
return t->nrec;
}

/*:46*//*47:*/
#line 1417 "diag180.w"

static void add_list(table_t*t,size_t ti,const table_t*src,size_t si,
int add,int limit)
{
uint32_t prev= 0,r= t->head[ti];
for(uint32_t s= src->head[si];s&&src->rec[s-1].wt<=limit;
s= src->rec[s-1].next){
int wt= src->rec[s-1].wt+add;
while(r&&t->rec[r-1].wt<wt){
prev= r;
r= t->rec[r-1].next;
}
if(r&&t->rec[r-1].wt==wt)
t->rec[r-1].c= addm(t->rec[r-1].c,src->rec[s-1].c);
else{
uint32_t q= new_record(t,wt,src->rec[s-1].c,r);
if(prev)
t->rec[prev-1].next= q;
else
t->head[ti]= q;
prev= q;
}
}
}

/*:47*//*49:*/
#line 1457 "diag180.w"

static int completion_possible(const state_t*s)
{
for(int j= 1;j<=s->num_classes;j++)
if(s->partner[j]!=0)
return 1;
return s->num_classes==0;
}

/*:49*//*57:*/
#line 1656 "diag180.w"

static inline int frow(int j)
{
return((umin+j-CUR_V)%2+2)%2?CUR_V-1:CUR_V;
}

/*:57*//*58:*/
#line 1662 "diag180.w"

static inline void reach(int h)
{
H0= max(H0,h);
}

static inline void sreach(int sc)
{
S0= min(S0,sc);
}

/*:58*//*59:*/
#line 1675 "diag180.w"

static inline int both_cost(int h,int j)
{
return max(VSTAR-h,min(j,W-1-j));
}

/*:59*//*60:*/
#line 1686 "diag180.w"

static int gap_score(int a,int ra,int b,int rb)
{
int best= 1<<20;
for(int j= a+1;j<=b;j++){
int t1= ra+1+(j-a-1),t2= rb+(b-j);
int c= both_cost(min(t1,t2),j);
best= min(best,c);
}
return best;
}

/*:60*//*61:*/
#line 1704 "diag180.w"

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
p1= max(p1,pk[t]);
s1= min(s1,sc[t]);
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

/*:61*//*62:*/
#line 1745 "diag180.w"

static inline int finish(int h,int s)
{
int top= max(VSTAR-h,0);
if(CUR_TOUCHES)
#line 33 "production.ch"
return top;
#line 1753 "diag180.w"
int sep= DEDGE+top;
return min(s,sep);
}

/*:62*//*63:*/
#line 1759 "diag180.w"

static int height_extra(int h,int s)
{
int ex= finish(h,s);
for(int o= 0;o<NOPT&&ex> 0;o++){
int e2= OPTD[o]+finish(max(OPTP[o],h),min(OPTS[o],s));
ex= min(ex,e2);
}
return ex;
}

/*:63*//*65:*/
#line 1988 "diag180.w"

static int cactus_bound(const state_t*s,int kink)
{
const uint8_t*label= s->label;
int touches= s->touches;
/*66:*/
#line 2004 "diag180.w"

int rs[MAXW],re[MAXW],rl[MAXW],R= 0;
int lastrun[MAXW+2];
for(int i= 0;i<W;i++)
if(label[i]){
if(i==0||!label[i-1]){
rs[R]= i;
rl[R]= label[i];
R++;
}
re[R-1]= i;
lastrun[label[i]]= R-1;
}
if(!R)
return 0;

/*:66*/
#line 1993 "diag180.w"

/*67:*/
#line 2020 "diag180.w"

H0= -1000;
S0= 1<<20;
NOPT= 0;
NCL= 0;
CUR_TOUCHES= touches;
for(int i= 0;i<W;i++)
if(label[i]){
reach(frow(i));
sreach(both_cost(frow(i),i));
}

/*:67*/
#line 1994 "diag180.w"

int cost= 0;
/*68:*/
#line 2039 "diag180.w"


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

/*:68*/
#line 1996 "diag180.w"

/*69:*/
#line 2078 "diag180.w"


int cr[MAXW],m= 0;
for(int t= 0;t<nt;t++)
if(s->partner[toplab[t]]!=0)
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

/*:69*/
#line 1997 "diag180.w"

/*70:*/
#line 2131 "diag180.w"

DEDGE= 0;
if(!touches){
int d1= rs[0],d2= W-1-re[R-1];
DEDGE= min(d1,d2);
}


if(NCL){
int ex= 1<<20;
for(int c= 0;c<NCL;c++){
int e2= CLD[c]+height_extra(max(CLP[c],H0),min(CLS[c],S0));
ex= min(ex,e2);
}
cost+= ex;
}else
cost+= height_extra(H0,S0);

#line 2217 "diag180.w"
/*:70*/
#line 1998 "diag180.w"

return cost;
}

/*:65*//*71:*/
#line 2224 "diag180.w"

static int bound(const state_t*s,int k,int kmin,int i,int last,int first_twin)
{
const uint8_t*label= s->label;
int lowerbound= max(2*(kmin-k),0);
if(!s->touches){
int dist= W;
for(int j= 0;j<W;j++)
if(label[j]){
int dd= min(j,W-1-j);
dist= min(dist,dd);
}
lowerbound= max(lowerbound,2*dist);
}
if(!first_twin||last){
if(!completion_possible(s))
return-1;
int kink= (!last&&!label[i]&&i+1<W&&label[i+1])?i+1:-1;
CUR_V= V0+k-1;
VSTAR= V0+kmin-1;
#line 111 "production.ch"
int eb= 2*cactus_bound(s,kink);
#line 2245 "diag180.w"
lowerbound= max(lowerbound,eb);
}
return lowerbound;
}

/*:71*//*74:*/
#line 2270 "diag180.w"

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

/*:74*//*75:*/
#line 2309 "diag180.w"

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

/*:75*//*77:*/
#line 2341 "diag180.w"

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

/*:77*//*78:*/
#line 2365 "diag180.w"

static void count_case(void)
{
/*79:*/
#line 2388 "diag180.w"

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
KW= (3*W+MAXC+1+63)/64;
NMAX_plus_1= NMAX+1;

int kmin= 0;
while(vspan(kmin)<W)
kmin++;

/*:79*/
#line 2368 "diag180.w"

/*80:*/
#line 2418 "diag180.w"

table_t A,B,*cur= &A,*nxt= &B;
tab_init(&A,1024);
tab_init(&B,1024);
tab_clear(cur);
state_t key;
int v0;

/*:80*/
#line 2369 "diag180.w"

/*81:*/
#line 2433 "diag180.w"

size_t t_init;
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
par[a+1]= f[j];
}
int touches= isc[0]||isc[W-1];
canon(f,par,touches,&key);
t_init= tab_get(cur,&key);
cur->head[t_init]= new_record(cur,nc,1,0);
}
v0= vm+1;
V0= v0;
}else{
memset(&key,0,sizeof(key));
t_init= tab_get(cur,&key);
cur->head[t_init]= new_record(cur,0,1,0);
v0= (M_v+1)/2;
V0= v0;
}
tab_finish(cur);

/*:81*/
#line 2370 "diag180.w"

int cells_v[MAXW];
for(int k= 1;cur->n;k++){
/*82:*/
#line 2481 "diag180.w"

int v= v0+k-1,nu= 0;
for(int i= 0;i<W;i++)
if(((umin+i-v)%2+2)%2==0)
cells_v[nu++]= i;
int first_twin= (E==1&&k==1);

/*:82*/
#line 2373 "diag180.w"

for(int ci= 0;ci<nu;ci++){
/*83:*/
#line 2488 "diag180.w"

int i= cells_v[ci];
int last= ci==nu-1;
tab_clear(nxt);
for(size_t si= 0;si<cur->n;si++){
/*84:*/
#line 2499 "diag180.w"

state_t sk;
unpack(cur->keys+si*KW,&sk);
for(int occ= 0;occ<=1;occ++){
if(occ&&k> 1&&sk.num_classes==0)
continue;
int status= step(&sk,i,occ,last,first_twin,&key);
if(status==DEAD)
continue;
if(status==COMPLETE){
record(k-1,cur,si,sk.touches);
continue;
}
if(key.num_classes==0&&!(first_twin&&!last))
continue;
int add= occ?2:0;
/*85:*/
#line 2529 "diag180.w"

uint64_t w[12];
pack_key(&key,w);
size_t ti= tab_find(nxt,w);
int lowerbound= 0;
if(ti!=NO_STATE)
lowerbound= nxt->lowerbound[ti];
else if(key.num_classes> 0){
lowerbound= bound(&key,k,kmin,i,last,first_twin);
if(lowerbound<0)
continue;
/*86:*/
#line 122 "production.ch"

if(!last&&!first_twin&&!key.label[i]&&i+1<W&&key.label[i+1]){
#line 2563 "diag180.w"

int i2= cells_v[ci+1],last2= ci+1==nu-1;
int la= 1<<20;
state_t k2;
for(int o2= 0;o2<=1&&la;o2++){
int st2= step(&key,i2,o2,last2,0,&k2);
if(st2==DEAD)
continue;

if(st2==COMPLETE&&!(key.touches&&vspan(k-1)>=W))
continue;
int val= 2*o2;
if(st2==OK){
if(k2.num_classes==0)
continue;
int lowerbound2= bound(&k2,k,kmin,i2,last2,0);
if(lowerbound2<0)
continue;
val+= lowerbound2;
}
la= min(la,val);
}
lowerbound= max(lowerbound,la);
}

/*:86*/
#line 2540 "diag180.w"

}
int limit= NMAX-lowerbound-add;
if(cur->rec[cur->head[si]-1].wt> limit)
continue;
if(ti==NO_STATE){
ti= tab_create(nxt,w);
nxt->lowerbound[ti]= lowerbound;
}
add_list(nxt,ti,cur,si,add,limit);

/*:85*/
#line 2515 "diag180.w"

}

/*:84*/
#line 2493 "diag180.w"

}
/*87:*/
#line 2588 "diag180.w"

table_t*x= cur;
cur= nxt;
nxt= x;
MAXSTATES= max(MAXSTATES,cur->n);
size_t by= tab_bytes(cur)+tab_bytes(nxt);
if(by> MAXBYTES)
MAXBYTES= by,MAXRECS= cur->nrec,MAXKW= KW;
tab_finish(cur);
if(!last&&(VERB>=2||(VERB==1&&wall()-T_LAST>=PSEC)))
stat_line("#cell",k,v,ci+1,nu,cur,nxt);

/*:87*/
#line 2495 "diag180.w"


/*:83*/
#line 2375 "diag180.w"

}
/*88:*/
#line 2608 "diag180.w"

if(REFLECT&&M_u==0&&cur->n){

tab_clear(nxt);
state_t sk,mirror;
for(size_t si= 0;si<cur->n;si++){
unpack(cur->keys+si*KW,&sk);
uint8_t f[MAXW];
for(int a= 0;a<W;a++)
f[a]= sk.label[W-1-a];
canon(f,sk.partner,sk.touches,&mirror);
const state_t*use= state_cmp(&mirror,&sk)<0?&mirror:&sk;
add_list(nxt,tab_get(nxt,use),cur,si,0,NMAX);
}
table_t*x= cur;
cur= nxt;
nxt= x;
tab_finish(cur);
}

/*:88*/
#line 2377 "diag180.w"

if(VERB>=1)
stat_line("#row ",k,v,nu,nu,cur,nxt);
if(vspan(k-1)>=W)
rect_line(k-1);
}
tab_free(&A);
tab_free(&B);
}

/*:78*/
#line 171 "diag180.w"

/*90:*/
#line 2648 "diag180.w"

int main(int argc,char**argv)
{
/*91:*/
#line 2666 "diag180.w"

int a= 1;
while(a<argc&&argv[a][0]=='-'){
if(!strcmp(argv[a],"-m")&&a+1<argc)
MOD= strtoull(argv[++a],0,10),a++;
else if(!strcmp(argv[a],"-o")&&a+1<argc)
RESFILE= argv[++a],a++;
#line 2681 "diag180.w"
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
/*73:*/
#line 2263 "diag180.w"

PORTRAIT= calloc(NMAX_plus_1,sizeof(cnt_t));
SQUARE= calloc(NMAX_plus_1,sizeof(cnt_t));
RECT= calloc(NMAX_plus_1,sizeof(cnt_t));

/*:73*/
#line 2706 "diag180.w"


/*:91*/
#line 2651 "diag180.w"

clock_t t0= clock();
T_START= T_LAST= wall();
/*93:*/
#line 2735 "diag180.w"

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
#line 161 "production.ch"
fprintf(f,"# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s; pack check %s\n",
REFLECT?"on":"off (-R)",CHECK?"on (-C)":"off");
#line 2752 "diag180.w"
fprintf(f,"# output: results file %s, progress -v %d -p %g s\n",
RESFILE,VERB,PSEC);
fprintf(f,"# sizes: counter %zu bytes, count record %zu bytes, key word 64 bits\n",
sizeof(cnt_t),sizeof(rec_t));
fflush(f);
}

/*:93*/
#line 2654 "diag180.w"

if(argc-a==4){
/*94:*/
#line 2761 "diag180.w"

M_x= atoi(argv[a+1]),M_y= atoi(argv[a+2]),umax= atoi(argv[a+3]);
count_case();
if(MOD)
printf("# modulus %llu; n #(uspan<vspan) #(uspan=vspan)\n",(unsigned long long)MOD);
else
printf("# modulus 2**"CNT_BITS"; n #(uspan<vspan) #(uspan=vspan)\n");
for(int n= 1;n<=NMAX;n++)
printf("%d %s %s\n",n,cs(PORTRAIT[n]),cs(SQUARE[n]));

/*:94*/
#line 2656 "diag180.w"

}else{
/*95:*/
#line 2774 "diag180.w"


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

/*:95*/
#line 2658 "diag180.w"

/*96:*/
#line 2804 "diag180.w"



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

/*:96*/
#line 2659 "diag180.w"

/*97:*/
#line 2828 "diag180.w"

printf("# n fix(r180) R180C(center in a cell) R180M(edge midpoint) R180V(vertex)\n");
for(int n= 1;n<=NMAX;n++){
cnt_t tot= addm(addm(out[n],out[NMAX_plus_1+n]),out[2*NMAX_plus_1+n]);
printf("%d %s %s %s %s\n",n,cs(tot),cs(out[n]),cs(out[NMAX_plus_1+n]),
cs(out[2*NMAX_plus_1+n]));
}

/*:97*/
#line 2660 "diag180.w"

}
/*98:*/
#line 2836 "diag180.w"

fprintf(stderr,"# peak memory %.1f MB (both tables; %zu count records, key %zu 64-bit words)\n",
MAXBYTES/1e6,MAXRECS,MAXKW);
fprintf(stderr,"# max states %zu, %.2fs\n",MAXSTATES,
(double)(clock()-t0)/CLOCKS_PER_SEC);

/*:98*/
#line 2662 "diag180.w"

return 0;
}

/*:90*/
#line 172 "diag180.w"


/*:2*/
