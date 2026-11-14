#define max(a,b) ((a) > (b) ?(a) :(b) ) 
#define min(a,b) ((a) <(b) ?(a) :(b) )  \
 \

#define MAXW 160
#define MAXLAB (MAXW+2) 
#define MAX_BYTES_PER_KEY 41
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
#line 651 "diag180.w"

typedef struct{
uint8_t label[MAXW];
uint8_t partner[MAXLAB];
uint8_t num_components;
uint8_t touches;
}state_t;

/*:13*//*19:*/
#line 781 "diag180.w"

enum{OK,DEAD,COMPLETE};

/*:19*//*24:*/
#line 903 "diag180.w"

enum{C_EMPTY,C_CONTINUE,C_SINGLETON,C_FIRST,C_MIDDLE,C_LAST};

/*:24*//*27:*/
#line 956 "diag180.w"

typedef struct{
uint8_t code[MAXW];
uint8_t crossing[MAXLAB];
uint8_t num_components;
uint8_t touches;
}coded_t;

/*:27*//*33:*/
#line 1107 "diag180.w"

#ifdef CNT128
typedef unsigned __int128 counter_t;
#define COUNTER_BITS "128"
#else
#line 1112 "diag180.w"
 typedef uint64_t counter_t;
#define COUNTER_BITS "64"
#endif
#line 1115 "diag180.w"

/*:33*//*37:*/
#line 1177 "diag180.w"

typedef struct{
counter_t count;
uint32_t next;
uint8_t weight;
}record_t;

typedef struct{
size_t num_states,state_capacity,hash_mask;
int decided_position;
uint8_t*keys;
uint32_t*first_record;
record_t*records;
size_t num_records,record_capacity;
uint32_t*hash;
uint8_t*lowerbound;
}table_t;

/*:37*/
#line 169 "diag180.w"

/*8:*/
#line 363 "diag180.w"

static int NMAX;
static int M_x,M_y;
static int M_u,M_v;
static int umin,umax;
static int W;
static int TWIN;
static int NMAX_plus_1;

/*:8*//*14:*/
#line 659 "diag180.w"

static int MAX_COMPONENTS;

/*:14*//*23:*/
#line 900 "diag180.w"

static int BYTES_PER_KEY;

/*:23*//*34:*/
#line 1116 "diag180.w"

static uint64_t MODULUS= 0;

/*:34*//*40:*/
#line 1275 "diag180.w"

static uint32_t HASH_BITS[MAX_BYTES_PER_KEY][256];

/*:40*//*58:*/
#line 1698 "diag180.w"

static int CURRENT_ROW,TARGET_ROW,FIRST_ROW;
#line 1701 "diag180.w"
static int FREE_HEIGHT,FREE_SCORE;
static int NUM_OPTIONS,OPTION_COST[2*MAXW+4],OPTION_PEAK[2*MAXW+4],OPTION_SCORE[2*MAXW+4];
static int NUM_MIDDLE,MIDDLE_GAP[MAXW],MIDDLE_PEAK[MAXW],MIDDLE_SCORE[MAXW];
static int SIDE_DISTANCE;
static int CURRENT_TOUCHES;
static int CURRENT_POSITION;

/*:58*//*74:*/
#line 2371 "diag180.w"

static counter_t*ROW_RESULTS;
static counter_t*PORTRAIT,*SQUARE;
static size_t MAX_STATES,MAX_BYTES,MAX_RECORDS,MAX_KEY_BYTES;

/*:74*//*78:*/
#line 2451 "diag180.w"

static int VERBOSITY= 1;
static double PROGRESS_SECONDS= 60;
static double START_TIME,LAST_REPORT_TIME;

/*:78*//*97:*/
#line 2885 "diag180.w"

static const uint8_t MIRRORED_CODE[6]= {C_EMPTY,C_CONTINUE,C_SINGLETON,C_LAST,C_MIDDLE,C_FIRST};

/*:97*//*100:*/
#line 2930 "diag180.w"

static int MERGE_MIRROR_STATES= 1;
#line 2933 "diag180.w"
static const char*RESULT_FILE_NAME= "results180.py";
#line 2935 "diag180.w"
static int CHECK_PACKING= 0;

/*:100*/
#line 170 "diag180.w"

/*9:*/
#line 391 "diag180.w"

static int vspan(int k){return 2*k+1-TWIN;}

/*:9*//*15:*/
#line 664 "diag180.w"

static int compare_states(const state_t*a,const state_t*b)
{

int c= memcmp(a->label,b->label,W);
if(c)
return c;
c= memcmp(a->partner+1,b->partner+1,a->num_components);
if(c)
return c;
return a->touches-b->touches;
}

/*:15*//*16:*/
#line 688 "diag180.w"

static inline int is_central(const uint8_t*partner,int a){return partner[a]==a;}
static inline int is_floating(const uint8_t*partner,int a){return partner[a]==0;}
static inline int is_paired(const uint8_t*partner,int a)
{return partner[a]!=0&&partner[a]!=a;}

static inline void relabel(uint8_t*label,int from,int to)
{
for(int j= 0;j<W;j++)
if(label[j]==from)
label[j]= to;
}

/*:16*//*17:*/
#line 715 "diag180.w"

static void join(uint8_t*label,uint8_t*partner,int a,int b)
{



if(a==b)
return;
int partner_a= partner[a],partner_b= partner[b];
relabel(label,b,a);
if(partner_a==b)
partner[a]= a;

else if(partner_a==a||partner_b==b){
if(partner_a!=0&&partner_a!=a)
relabel(label,partner_a,a);
if(partner_b!=0&&partner_b!=b)
relabel(label,partner_b,a);
partner[a]= a;
}else if(partner_a!=0){
if(partner_b!=0)
relabel(label,partner_b,partner_a);
}else if(partner_b!=0){
partner[a]= partner_b;
partner[partner_b]= a;
}
}

/*:17*//*18:*/
#line 745 "diag180.w"

static void canonical_state(const uint8_t*label,const uint8_t*partner,int touches,state_t*result)
{


uint8_t new_label[MAXLAB],old_label[MAXLAB];
int num_components= 0;
memset(new_label,0,sizeof(new_label));
for(int i= 0;i<W;i++){
if(label[i]&&!new_label[label[i]]){
new_label[label[i]]= ++num_components;
old_label[num_components]= label[i];
}
result->label[i]= new_label[label[i]];
}
for(int j= 1;j<=num_components;j++)
result->partner[j]= new_label[partner[old_label[j]]];
result->num_components= num_components;
result->touches= touches;
}

/*:18*//*20:*/
#line 784 "diag180.w"

static int retire(uint8_t*label,uint8_t*partner,int i)
{

int retired_label= label[i];
if(!retired_label)
return OK;
label[i]= 0;
int others_remain= 0;
for(int j= 0;j<W;j++)
if(label[j]){
others_remain= 1;
if(label[j]==retired_label)
return OK;
}
if(is_central(partner,retired_label))
return others_remain?DEAD:COMPLETE;
if(is_floating(partner,retired_label))
return DEAD;
partner[partner[retired_label]]= 0;
return OK;
}

/*:20*//*21:*/
#line 822 "diag180.w"

static int step(const state_t*state,int i,int occupied,int last_in_row,state_t*result)
{



uint8_t label[MAXW],partner[MAXLAB];
memcpy(label,state->label,W);
memcpy(partner,state->partner,state->num_components+1);
if(occupied){
int new_label= state->num_components+1;
partner[new_label]= 0;
label[i]= new_label;
for(int neighbor= i-1;neighbor<=i+1;neighbor+= 2)
if(neighbor>=0&&neighbor<W&&label[neighbor])
join(label,partner,new_label,label[neighbor]);
}
int status= OK;
if(i-1>=0)
status= retire(label,partner,i-1);
if(status==OK&&last_in_row&&i+1<W)
status= retire(label,partner,i+1);
if(status!=OK)
return status;
canonical_state(label,partner,state->touches||(occupied&&(i==0||i==W-1)),result);
return OK;
}

/*:21*//*26:*/
#line 928 "diag180.w"

static inline void put_bits(uint8_t*bytes,int*bit_position,int value,int num_bits)
{

int byte= *bit_position>>3,offset= *bit_position&7;
bytes[byte]|= value<<offset;
if(offset+num_bits> 8)
bytes[byte+1]|= value>>(8-offset);
*bit_position+= num_bits;
}

static inline int get_bits(const uint8_t*bytes,int*bit_position,int num_bits)
{
int byte= *bit_position>>3,offset= *bit_position&7;
int value= bytes[byte]>>offset;
if(offset+num_bits> 8)
value|= bytes[byte+1]<<(8-offset);
*bit_position+= num_bits;
return value&((1<<num_bits)-1);
}

/*:26*//*28:*/
#line 966 "diag180.w"

static void encode(const state_t*s,coded_t*coded)
{

int first_run[MAXLAB],last_run[MAXLAB];
int num_components= 0;
for(int i= 0;i<W;i++)
if(s->label[i]&&(i==0||!s->label[i-1])){
int lab= s->label[i];
if(lab> num_components){
num_components= lab;
first_run[lab]= i;
}
last_run[lab]= i;
}
for(int i= 0;i<W;i++){
int lab= s->label[i];
if(!lab)
coded->code[i]= C_EMPTY;
else if(i> 0&&s->label[i-1])
coded->code[i]= C_CONTINUE;
else if(first_run[lab]==i)
coded->code[i]= last_run[lab]==i?C_SINGLETON:C_FIRST;
else
coded->code[i]= last_run[lab]==i?C_LAST:C_MIDDLE;
}
for(int j= 1;j<=num_components;j++)
coded->crossing[j]= s->partner[j]!=0;
coded->num_components= num_components;
coded->touches= s->touches;
}

/*:28*//*29:*/
#line 1002 "diag180.w"

static void decode(const coded_t*coded,state_t*s)
{

int stack[MAXLAB],stack_size= 0,num_components= 0;
for(int i= 0;i<W;i++){
switch(coded->code[i]){
case C_EMPTY:s->label[i]= 0;break;
case C_CONTINUE:s->label[i]= s->label[i-1];break;
case C_SINGLETON:s->label[i]= ++num_components;break;
case C_FIRST:s->label[i]= ++num_components;stack[stack_size++]= num_components;break;
case C_MIDDLE:s->label[i]= stack[stack_size-1];break;
case C_LAST:s->label[i]= stack[--stack_size];break;
}
}
stack_size= 0;
for(int j= 1;j<=num_components;j++)
if(coded->crossing[j])
stack[stack_size++]= j;
else
s->partner[j]= 0;
for(int a= 0,b= stack_size-1;a<=b;a++,b--){
s->partner[stack[a]]= stack[b];
s->partner[stack[b]]= stack[a];
}
s->num_components= num_components;
s->touches= coded->touches;
}

/*:29*//*30:*/
#line 1040 "diag180.w"

static inline int can_hold_frontier_cell(int j,int decided_position)
{
return(((j^decided_position)&1)!=0)==(j> decided_position);
}

/*:30*//*31:*/
#line 1049 "diag180.w"

static void compress(const coded_t*coded,int decided_position,uint8_t*bytes)
{

int pos= 0;
memset(bytes,0,BYTES_PER_KEY);
for(int i= 0;i<W;i++)
if(can_hold_frontier_cell(i,decided_position))
put_bits(bytes,&pos,coded->code[i],3);
for(int j= 1;j<=coded->num_components;j++)
put_bits(bytes,&pos,coded->crossing[j],1);
put_bits(bytes,&pos,coded->touches,1);
}

static void expand(const uint8_t*bytes,int decided_position,coded_t*coded)
{
int pos= 0,m= 0;
for(int i= 0;i<W;i++){
coded->code[i]= can_hold_frontier_cell(i,decided_position)?get_bits(bytes,&pos,3):C_EMPTY;
if(coded->code[i]==C_SINGLETON||coded->code[i]==C_FIRST)
m++;
}
for(int j= 1;j<=m;j++)
coded->crossing[j]= get_bits(bytes,&pos,1);
coded->num_components= m;
coded->touches= get_bits(bytes,&pos,1);
}

/*:31*//*32:*/
#line 1077 "diag180.w"

static void pack(const state_t*s,int decided_position,uint8_t*bytes)
{
coded_t coded;
encode(s,&coded);
compress(&coded,decided_position,bytes);
}

static void unpack(const uint8_t*bytes,int decided_position,state_t*s)
{
coded_t coded;
expand(bytes,decided_position,&coded);
decode(&coded,s);
}

/*:32*//*35:*/
#line 1119 "diag180.w"

static inline counter_t add_mod(counter_t a,counter_t b)
{

if(!MODULUS)
return a+b;
counter_t s= a+b;
if(s<a||s>=MODULUS)
s-= MODULUS;
return s;
}

/*:35*//*36:*/
#line 1135 "diag180.w"

static const char*counter_to_string(counter_t x)
{
static char buffer[8][44];
static int current_buffer;
char*end= buffer[current_buffer= (current_buffer+1)&7]+43,*digit= end;
*digit= 0;
do
*--digit= '0'+(int)(x%10),x/= 10;
while(x);
return digit;
}

/*:36*//*38:*/
#line 1195 "diag180.w"

static void*checked_realloc(void*p,size_t s)
{
p= realloc(p,s);
if(!p){
fprintf(stderr,"out of memory\n");
exit(1);
}
return p;
}

/*:38*//*39:*/
#line 1206 "diag180.w"

static void table_allocate(table_t*t,size_t capacity)
{

t->num_states= 0;
t->state_capacity= capacity;
size_t hash_size= 16;
while(hash_size<2*capacity)
hash_size<<= 1;
t->hash_mask= hash_size-1;
t->keys= checked_realloc(0,capacity*BYTES_PER_KEY);
t->first_record= checked_realloc(0,capacity*4);
t->hash= 0;
t->lowerbound= 0;
t->record_capacity= 2*capacity;
t->num_records= 0;
t->records= checked_realloc(0,t->record_capacity*sizeof(record_t));
}

static void table_free(table_t*t)
{
free(t->keys);
free(t->first_record);
free(t->records);
free(t->hash);
free(t->lowerbound);
}

static size_t table_bytes(const table_t*t)
{
return t->num_states*(BYTES_PER_KEY+4)+t->num_records*sizeof(record_t)
+(t->hash?(t->hash_mask+1)*4+t->num_states:0);
}

/*:39*//*41:*/
#line 1278 "diag180.w"

static inline uint32_t hash_of_key(const uint8_t*key)
{
uint32_t h= HASH_BITS[0][key[0]];
for(int l= 1;l<BYTES_PER_KEY;l++)
h+= HASH_BITS[l][key[l]];
return h;
}

/*:41*//*44:*/
#line 1338 "diag180.w"

static void table_initialize(table_t*t,int decided_position)
{

t->decided_position= decided_position;
t->num_states= 0;
t->num_records= 0;
t->hash= calloc(t->hash_mask+1,sizeof(uint32_t));
if(!t->hash){
fprintf(stderr,"out of memory\n");
exit(1);
}
t->lowerbound= checked_realloc(0,t->state_capacity);
}

static void table_finish(table_t*t)
{
free(t->hash);
t->hash= 0;
free(t->lowerbound);
t->lowerbound= 0;
}

/*:44*//*45:*/
#line 1369 "diag180.w"

static size_t HASH_SLOT;

static size_t table_find(const table_t*t,const uint8_t*packed_key)
{

size_t slot= hash_of_key(packed_key)&t->hash_mask;
while(t->hash[slot]){
size_t i= t->hash[slot]-1;
if(!memcmp(t->keys+i*BYTES_PER_KEY,packed_key,BYTES_PER_KEY))
return i;
slot= (slot+1)&t->hash_mask;
}
HASH_SLOT= slot;
return NO_STATE;
}

static size_t table_create(table_t*t,const uint8_t*packed_key)
{

size_t slot= HASH_SLOT;
if(t->num_states==t->state_capacity){
/*43:*/
#line 1307 "diag180.w"

size_t capacity= t->state_capacity*2;
t->keys= checked_realloc(t->keys,capacity*BYTES_PER_KEY);
t->first_record= checked_realloc(t->first_record,capacity*4);
t->lowerbound= checked_realloc(t->lowerbound,capacity);
size_t hash_size= (t->hash_mask+1)*2;
free(t->hash);
t->hash= calloc(hash_size,sizeof(uint32_t));
if(!t->hash){
fprintf(stderr,"out of memory\n");
exit(1);
}
t->hash_mask= hash_size-1;
t->state_capacity= capacity;
for(size_t state= 0;state<t->num_states;state++){
size_t state_slot= hash_of_key(t->keys+state*BYTES_PER_KEY)&t->hash_mask;
while(t->hash[state_slot])
state_slot= (state_slot+1)&t->hash_mask;
t->hash[state_slot]= state+1;
}

/*:43*/
#line 1391 "diag180.w"

slot= hash_of_key(packed_key)&t->hash_mask;
while(t->hash[slot])
slot= (slot+1)&t->hash_mask;
}
size_t i= t->num_states++;
t->hash[slot]= i+1;
memcpy(t->keys+i*BYTES_PER_KEY,packed_key,BYTES_PER_KEY);
t->first_record[i]= 0;
t->lowerbound[i]= 0;
return i;
}

/*:45*//*46:*/
#line 1406 "diag180.w"

static void pack_and_check(const state_t*s,int decided_position,uint8_t*packed_key)
{
pack(s,decided_position,packed_key);
if(CHECK_PACKING){
state_t unpacked;
unpack(packed_key,decided_position,&unpacked);
if(unpacked.num_components!=s->num_components||compare_states(&unpacked,s)){
fprintf(stderr,"pack/unpack mismatch\n");
exit(2);
}
}
}

/*:46*//*47:*/
#line 1422 "diag180.w"

static size_t table_find_or_create(table_t*t,const state_t*s)
{
uint8_t packed_key[MAX_BYTES_PER_KEY];
pack_and_check(s,t->decided_position,packed_key);
size_t i= table_find(t,packed_key);
return i!=NO_STATE?i:table_create(t,packed_key);
}

/*:47*//*48:*/
#line 1435 "diag180.w"

static uint32_t new_record(table_t*t,int weight,counter_t count,uint32_t next)
{

if(t->num_records==t->record_capacity){
if(t->record_capacity>=0xFFFFFFFFu){
fprintf(stderr,"too many count records for 32-bit links\n");
exit(1);
}
t->record_capacity= min(2*t->record_capacity,0xFFFFFFFFu);
t->records= checked_realloc(t->records,t->record_capacity*sizeof(record_t));
}
record_t*q= &t->records[t->num_records++];
q->count= count;
q->weight= weight;
q->next= next;
return t->num_records;
}

/*:48*//*49:*/
#line 1465 "diag180.w"

static void add_list(table_t*target,size_t target_state,const table_t*source,size_t source_state,
int added_weight,int weight_limit)
{


uint32_t previous= 0,t_index= target->first_record[target_state];
for(uint32_t s_index= source->first_record[source_state];s_index&&source->records[s_index-1].weight<=weight_limit;
s_index= source->records[s_index-1].next){
int weight= source->records[s_index-1].weight+added_weight;
while(t_index&&target->records[t_index-1].weight<weight){
previous= t_index;
t_index= target->records[t_index-1].next;
}
if(t_index&&target->records[t_index-1].weight==weight)
target->records[t_index-1].count= add_mod(target->records[t_index-1].count,source->records[s_index-1].count);
else{
uint32_t new_one= new_record(target,weight,source->records[s_index-1].count,t_index);
if(previous)
target->records[previous-1].next= new_one;
else
target->first_record[target_state]= new_one;
previous= new_one;
}
}
}

/*:49*//*59:*/
#line 1712 "diag180.w"

static inline int frontier_row(int j)
{


return j<=CURRENT_POSITION?CURRENT_ROW:CURRENT_ROW-1;
}

/*:59*//*60:*/
#line 1720 "diag180.w"

static inline void reach_height(int h)
{
FREE_HEIGHT= max(FREE_HEIGHT,h);
}

static inline void reach_score(int sc)
{
FREE_SCORE= min(FREE_SCORE,sc);
}

/*:60*//*61:*/
#line 1733 "diag180.w"

static inline int side_and_top_cost(int h,int j)
{
return max(TARGET_ROW-h,min(j,W-1-j));
}

/*:61*//*62:*/
#line 1758 "diag180.w"

static int gap_score(int a,int row_a,int b,int row_b)
{

int best= 1<<20;
for(int j= a+1;j<b;j++){
int row_from_left= row_a+(j-a),row_from_right= row_b+(b-j);
int c= side_and_top_cost(min(row_from_left,row_from_right),j);
best= min(best,c);
}
return best;
}

/*:62*//*63:*/
#line 1777 "diag180.w"

static int cycle_cost(const int*gap,const int*peak,const int*score,int num_gaps,int*max1_out,int*peak_out,int*score_out)
{

int sum= 0,max1= -1,num_max= 0,max2= 0,peak_of_max= -1000,score_of_max= 1<<20;
for(int t= 0;t<num_gaps;t++){
sum+= gap[t];
if(gap[t]> max1){
if(max1>=0)
max2= max1;
max1= gap[t];
num_max= 1;
peak_of_max= peak[t];
score_of_max= score[t];
}else if(gap[t]==max1){
num_max++;
peak_of_max= max(peak_of_max,peak[t]);
score_of_max= min(score_of_max,score[t]);
}else if(gap[t]> max2)
max2= gap[t];
}
for(int t= 0;t<num_gaps;t++)
if(num_max>=2||gap[t]!=max1){
reach_height(peak[t]);
reach_score(score[t]);
}
if(num_max==1){
OPTION_COST[NUM_OPTIONS]= max1-max2;
OPTION_SCORE[NUM_OPTIONS]= score_of_max;
OPTION_PEAK[NUM_OPTIONS++]= peak_of_max;
}
if(max1_out){
*max1_out= max1;
*peak_out= peak_of_max;
*score_out= score_of_max;
}
return sum-max1;
}

/*:63*//*64:*/
#line 1819 "diag180.w"

static inline int finish_cost(int height,int score)
{

int missing_rows= max(TARGET_ROW-height,0);
if(CURRENT_TOUCHES)
#line 33 "production.ch"
return missing_rows;
#line 1828 "diag180.w"
int separately= SIDE_DISTANCE+missing_rows;
return min(score,separately);
}

/*:64*//*65:*/
#line 1834 "diag180.w"

static int height_term(int height,int score)
{
int extra= finish_cost(height,score);
for(int option= 0;option<NUM_OPTIONS&&extra> 0;option++){
int extra_with_option= OPTION_COST[option]+finish_cost(max(OPTION_PEAK[option],height),min(OPTION_SCORE[option],score));
extra= min(extra,extra_with_option);
}
return extra;
}

/*:65*//*67:*/
#line 2098 "diag180.w"

static int cactus_bound(const state_t*state,int kink)
{


const uint8_t*label= state->label;
int touches= state->touches;
/*68:*/
#line 2116 "diag180.w"

int run_start[MAXW],run_end[MAXW],run_label[MAXW],num_runs= 0;
int last_run[MAXW+2];
for(int i= 0;i<W;i++)
if(label[i]){
if(i==0||!label[i-1]){
run_start[num_runs]= i;
run_label[num_runs]= label[i];
num_runs++;
}
run_end[num_runs-1]= i;
last_run[label[i]]= num_runs-1;
}
if(!num_runs)
return 0;

/*:68*/
#line 2105 "diag180.w"

/*69:*/
#line 2132 "diag180.w"

FREE_HEIGHT= -1000;
FREE_SCORE= 1<<20;
NUM_OPTIONS= 0;
NUM_MIDDLE= 0;
CURRENT_TOUCHES= touches;
for(int i= 0;i<W;i++)
if(label[i]){
reach_height(frontier_row(i));
reach_score(side_and_top_cost(frontier_row(i),i));
}

/*:69*/
#line 2106 "diag180.w"

int cost= 0;
/*70:*/
#line 2151 "diag180.w"


int pocket_label[MAXW],pocket_first_gap[MAXW],num_open_pockets= 0;
int pocket_gap[MAXW],pocket_peak[MAXW],pocket_score[MAXW],num_pocket_gaps= 0;
int top_label[MAXW],top_gap[MAXW],top_peak[MAXW],top_score[MAXW],num_top= 0;
for(int j= 0;j<num_runs;j++){
int component= run_label[j];
if(j> 0){
int a= run_end[j-1],b= run_start[j];
int gap= b-a-1+(b==kink);

int peak= (gap-1+(b==kink)+frontier_row(a)+frontier_row(b)+2)/2;
int score= touches?1<<20:gap_score(a,frontier_row(a),b==kink?b+2:b,frontier_row(b));
if(num_open_pockets){
pocket_gap[num_pocket_gaps]= gap;
pocket_score[num_pocket_gaps]= score;
pocket_peak[num_pocket_gaps++]= peak;
}else{
top_gap[num_top-1]= gap;
top_score[num_top-1]= score;
top_peak[num_top-1]= peak;
}
}
if(num_open_pockets&&pocket_label[num_open_pockets-1]==component){
num_open_pockets--;
cost+= cycle_cost(pocket_gap+
pocket_first_gap[num_open_pockets],pocket_peak+
pocket_first_gap[num_open_pockets],pocket_score+
pocket_first_gap[num_open_pockets],
num_pocket_gaps-pocket_first_gap[num_open_pockets],0,0,0);
num_pocket_gaps= pocket_first_gap[num_open_pockets];
}else if(!num_open_pockets)
top_label[num_top++]= component;
if(j<last_run[component]){
pocket_label[num_open_pockets]= component;
pocket_first_gap[num_open_pockets]= num_pocket_gaps;
num_open_pockets++;
}
}

/*:70*/
#line 2108 "diag180.w"

/*71:*/
#line 2194 "diag180.w"


int crossing[MAXW],num_crossing= 0;
for(int t= 0;t<num_top;t++)
if(state->partner[top_label[t]]!=0)
crossing[num_crossing++]= t;
if(num_crossing){
for(int t= 0;t<crossing[0];t++){
cost+= top_gap[t];
reach_height(top_peak[t]);
reach_score(top_score[t]);
}
for(int t= crossing[num_crossing-1];t<num_top-1;t++){
cost+= top_gap[t];
reach_height(top_peak[t]);
reach_score(top_score[t]);
}
int num_pairs= num_crossing/2;
NUM_MIDDLE= 0;
for(int i= 0;i<num_pairs;i++){
int gap[2*MAXW],peak[2*MAXW],score[2*MAXW],num_gaps= 0,max1,peak_of_max,score_of_max;
int a= crossing[i],b= (i<num_pairs-1)?crossing[i+1]:crossing[num_crossing-1-i];
for(int t= a;t<b;t++){
gap[num_gaps]= top_gap[t];
score[num_gaps]= top_score[t];
peak[num_gaps++]= top_peak[t];
}
if(i<num_pairs-1)
for(int t= crossing[num_crossing-2-i];t<crossing[num_crossing-1-i];t++){
gap[num_gaps]= top_gap[t];
score[num_gaps]= top_score[t];
peak[num_gaps++]= top_peak[t];
}
cost+= cycle_cost(gap,peak,score,num_gaps,&max1,&peak_of_max,&score_of_max);
if(i==num_pairs-1){



MIDDLE_GAP[NUM_MIDDLE]= max1;
MIDDLE_SCORE[NUM_MIDDLE]= score_of_max;
MIDDLE_PEAK[NUM_MIDDLE++]= peak_of_max;
}
}
if(num_crossing&1)
NUM_MIDDLE= 0;
}else{
for(int t= 0;t<num_top-1;t++){
cost+= top_gap[t];
reach_height(top_peak[t]);
reach_score(top_score[t]);
}
}

/*:71*/
#line 2109 "diag180.w"

/*72:*/
#line 2247 "diag180.w"

SIDE_DISTANCE= 0;
if(!touches){
int left_distance= run_start[0],right_distance= W-1-run_end[num_runs-1];
SIDE_DISTANCE= min(left_distance,right_distance);
}



if(NUM_MIDDLE){
int extra= 1<<20;
for(int c= 0;c<NUM_MIDDLE;c++){
int extra_this_way= MIDDLE_GAP[c]+height_term(max(MIDDLE_PEAK[c],FREE_HEIGHT),min(MIDDLE_SCORE[c],FREE_SCORE));
extra= min(extra,extra_this_way);
}
cost+= extra;
}else
cost+= height_term(FREE_HEIGHT,FREE_SCORE);

#line 2329 "diag180.w"
/*:72*/
#line 2110 "diag180.w"

return cost;
}

/*:67*//*73:*/
#line 2335 "diag180.w"

static int bound(const state_t*state,int k,int min_rows,int i,int last_in_row)
{


const uint8_t*label= state->label;
int lowerbound= max(2*(min_rows-k),0);
if(!state->touches){
int side_distance= W;
for(int j= 0;j<W;j++)
if(label[j]){
int distance= min(j,W-1-j);
side_distance= min(side_distance,distance);
}
lowerbound= max(lowerbound,2*side_distance);
}
/*51:*/
#line 1506 "diag180.w"

{
int some_component_is_crossing= 0;
for(int j= 1;j<=state->num_components;j++)
if(state->partner[j]!=0)
some_component_is_crossing= 1;
if(!some_component_is_crossing)
return-1;
}

/*:51*/
#line 2351 "diag180.w"

int kink= (!last_in_row&&!label[i]&&i+1<W&&label[i+1])?i+1:-1;
CURRENT_ROW= FIRST_ROW+k-1;
CURRENT_POSITION= i;
TARGET_ROW= FIRST_ROW+min_rows-1;
#line 110 "production.ch"
int connection_bound= 2*cactus_bound(state,kink);
#line 2361 "diag180.w"
lowerbound= max(lowerbound,connection_bound);
return lowerbound;
}

/*:73*//*77:*/
#line 2425 "diag180.w"

static void print_row_results(int top_row_number)
{
static const char*type_name[3]= {
"R180C center in a cell",
"R180M center on the midpoint of a vertical edge",
"R180V center at a vertex"};
int height= vspan(top_row_number);
int top_row= TWIN?(M_v+1)/2+top_row_number-1:M_v/2+top_row_number;
if(VERBOSITY>=1){
fprintf(stderr,"#rect %s (Mx=%d My=%d): width %d (u=%d..%d) x height %d (v=%d..%d):",
type_name[M_x+M_y],M_x,M_y,W,umin,umin+W-1,height,M_v-top_row,top_row);
for(int n= 0;n<=NMAX;n++)
if(ROW_RESULTS[n])
fprintf(stderr," n=%d:%s",n,counter_to_string(ROW_RESULTS[n]));
fprintf(stderr,"\n");
fflush(stderr);
}
}

/*:77*//*79:*/
#line 2456 "diag180.w"

static double wall_time(void)
{
struct timespec ts;
clock_gettime(CLOCK_MONOTONIC,&ts);
return ts.tv_sec+1e-9*ts.tv_nsec;
}
static void print_status_line(const char*tag,int k,int v,int cell_number,int num_cells,
const table_t*current,const table_t*next)
{
double now= wall_time();
LAST_REPORT_TIME= now;
counter_t done= 0;
for(int n= 0;n<=NMAX;n++)
done= add_mod(done,add_mod(add_mod(PORTRAIT[n],SQUARE[n]),ROW_RESULTS[n]));
fprintf(stderr,"%s Mx=%d My=%d umax=%d W=%d row k=%d v=%d cell %d/%d: "
"%zu states, %zu records, %.1f MB, done %s, %.1fs wall, %.1fs cpu\n",
tag,M_x,M_y,umax,W,k,v,cell_number,num_cells,current->num_states,current->num_records,
(table_bytes(current)+table_bytes(next))/1e6,counter_to_string(done),
now-START_TIME,(double)clock()/CLOCKS_PER_SEC);
fflush(stderr);
}

/*:79*//*80:*/
#line 2481 "diag180.w"

static void count_case(void)
{

/*82:*/
#line 2519 "diag180.w"

M_v= M_x+M_y;
M_u= M_x-M_y;
umin= M_u-umax;
W= umax-umin+1;
TWIN= M_v&1;
memset(PORTRAIT,0,NMAX_plus_1*sizeof(counter_t));
memset(SQUARE,0,NMAX_plus_1*sizeof(counter_t));
memset(ROW_RESULTS,0,NMAX_plus_1*sizeof(counter_t));
if(W==1){

if(M_x==0&&M_y==0){
SQUARE[1]= ROW_RESULTS[1]= 1;
print_row_results(0);
}
return;
}
if(W> MAXW){
fprintf(stderr,"strip too wide\n");
exit(1);
}
MAX_COMPONENTS= (W+1)/2+1;
BYTES_PER_KEY= (3*(W/2+1)+MAX_COMPONENTS+1+7)/8;
NMAX_plus_1= NMAX+1;

int min_rows= 0;
while(vspan(min_rows)<W)
min_rows++;

/*:82*/
#line 2485 "diag180.w"

/*83:*/
#line 2550 "diag180.w"

table_t A,B,*current= &A,*next= &B;
table_allocate(&A,1024);
table_allocate(&B,1024);
state_t new_state;
int first_row;

/*:83*/
#line 2486 "diag180.w"

/*84:*/
#line 2564 "diag180.w"

int choice[MAXW],num_choices= 0,max_chosen;
if(TWIN==0){
/*86:*/
#line 2621 "diag180.w"

int vm= M_v/2;
for(int i= 0;i<W;i++)
if(((umin+i-vm)%2+2)%2==0&&2*(umin+i)>=M_u)
choice[num_choices++]= i;
max_chosen= (NMAX-2*min_rows+1)/2;
first_row= vm+1;
FIRST_ROW= first_row;

/*:86*/
#line 2567 "diag180.w"

}else{
/*88:*/
#line 2656 "diag180.w"

first_row= (M_v+1)/2;
FIRST_ROW= first_row;
for(int i= 0;i<W;i++)
if(((umin+i-first_row)%2+2)%2==0)
choice[num_choices++]= i;
max_chosen= (NMAX-2*(min_rows-1))/2;

/*:88*/
#line 2569 "diag180.w"

}
uint8_t chosen[MAXW];
int num_chosen= 0,done= 0;
memset(chosen,0,num_choices);
table_initialize(current,choice[num_choices-1]);
while(!done){
/*85:*/
#line 2598 "diag180.w"

{
int j= 0;
if(num_chosen>=max_chosen)
while(j<num_choices&&!chosen[j])
j++;
while(j<num_choices&&chosen[j]){
chosen[j++]= 0;
num_chosen--;
}
if(j==num_choices)
done= 1;
else{
chosen[j]= 1;
num_chosen++;
}
}

/*:85*/
#line 2576 "diag180.w"

if(done)
break;
if(TWIN==0){
/*87:*/
#line 2633 "diag180.w"

uint8_t label[MAXW],partner[MAXLAB];
int num_cells= 0;
memset(label,0,W);
for(int b= 0;b<num_choices;b++)
if(chosen[b])
label[choice[b]]= label[W-1-choice[b]]= 1;
for(int i= 0;i<W;i++)
if(label[i])
label[i]= ++num_cells;
for(int i= 0;i<W;i++)
if(label[i])
partner[label[i]]= label[W-1-i];
if(num_cells+2*min_rows<=NMAX){
canonical_state(label,partner,label[0]!=0,&new_state);
size_t target_state= table_find_or_create(current,&new_state);
current->first_record[target_state]= new_record(current,num_cells,1,0);
}

/*:87*/
#line 2580 "diag180.w"

}else{
/*89:*/
#line 2683 "diag180.w"

uint8_t band[MAXW+1];
int is_twin[MAXW],run_label[MAXW],num_runs= 0,num_components= 0;
memset(band,0,W+1);
memset(is_twin,0,W*sizeof(int));
for(int b= 0;b<num_choices;b++){
is_twin[choice[b]]= 1;
if(chosen[b])
band[choice[b]]= band[W-1-choice[b]]= 1;
}
for(int q= 0;q<W;q++){
new_state.label[q]= 0;
if(!band[q])
continue;
if(q==0||!band[q-1])
run_label[num_runs++]= 0;
if(is_twin[q]){
if(!run_label[num_runs-1])
run_label[num_runs-1]= ++num_components;
new_state.label[q]= run_label[num_runs-1];
}
}
for(int r= 0;r<num_runs;r++)
if(run_label[r])
new_state.partner[run_label[r]]= run_label[num_runs-1-r];
new_state.num_components= num_components;
new_state.touches= band[0];
int lowerbound= bound(&new_state,1,min_rows,W-1,1);
if(lowerbound>=0&&2*num_chosen+lowerbound<=NMAX){
size_t target_state= table_find_or_create(current,&new_state);
current->first_record[target_state]= new_record(current,2*num_chosen,1,0);
}

/*:89*/
#line 2582 "diag180.w"

}
}
table_finish(current);

/*:84*/
#line 2487 "diag180.w"

int row_cells[MAXW];
for(int k= 1+TWIN;current->num_states;k++){
memset(ROW_RESULTS,0,NMAX_plus_1*sizeof(counter_t));
/*90:*/
#line 2719 "diag180.w"

int v= first_row+k-1,num_row_cells= 0;
for(int i= 0;i<W;i++)
if(((umin+i-v)%2+2)%2==0)
row_cells[num_row_cells++]= i;

/*:90*/
#line 2491 "diag180.w"

for(int cell_number= 0;cell_number<num_row_cells;cell_number++){
/*91:*/
#line 2725 "diag180.w"

int i= row_cells[cell_number];
int last_in_row= cell_number==num_row_cells-1;
table_initialize(next,i);
for(size_t source_state= 0;source_state<current->num_states;source_state++){
/*92:*/
#line 2736 "diag180.w"

state_t old_state;
unpack(current->keys+source_state*BYTES_PER_KEY,current->decided_position,&old_state);
for(int occupied= 0;occupied<=1;occupied++){
int status= step(&old_state,i,occupied,last_in_row,&new_state);
if(status==DEAD)
continue;
if(status==COMPLETE){
/*76:*/
#line 2391 "diag180.w"

if(old_state.touches&&vspan(k-1)>=W)
for(uint32_t r= current->first_record[source_state];r;r= current->records[r-1].next){
int weight= current->records[r-1].weight;
ROW_RESULTS[weight]= add_mod(ROW_RESULTS[weight],current->records[r-1].count);
}

/*:76*/
#line 2744 "diag180.w"

continue;
}
int added_weight= occupied?2:0;
/*93:*/
#line 2761 "diag180.w"

uint8_t packed_key[MAX_BYTES_PER_KEY];
pack_and_check(&new_state,next->decided_position,packed_key);
size_t target_state= table_find(next,packed_key);
int lowerbound;
if(target_state!=NO_STATE)
lowerbound= next->lowerbound[target_state];
else{
lowerbound= bound(&new_state,k,min_rows,i,last_in_row);
if(lowerbound<0)
continue;
/*94:*/
#line 121 "production.ch"

if(!last_in_row&&!new_state.label[i]&&i+1<W&&new_state.label[i+1]){
#line 2795 "diag180.w"

int next_i= row_cells[cell_number+1],next_is_last= cell_number+1==num_row_cells-1;
int lookahead_bound= 1<<20;
state_t state_ahead;
for(int next_occupied= 0;next_occupied<=1&&lookahead_bound;next_occupied++){
int status_ahead= step(&new_state,next_i,next_occupied,next_is_last,&state_ahead);
if(status_ahead==DEAD)
continue;

if(status_ahead==COMPLETE&&!(new_state.touches&&vspan(k-1)>=W))
continue;
int bound_this_way= 2*next_occupied;
if(status_ahead==OK){
int lowerbound_ahead= bound(&state_ahead,k,min_rows,next_i,next_is_last);
if(lowerbound_ahead<0)
continue;
bound_this_way+= lowerbound_ahead;
}
lookahead_bound= min(lookahead_bound,bound_this_way);
}
lowerbound= max(lowerbound,lookahead_bound);
}

/*:94*/
#line 2772 "diag180.w"

}
int limit= NMAX-lowerbound-added_weight;
if(current->records[current->first_record[source_state]-1].weight> limit)
continue;
if(target_state==NO_STATE){
target_state= table_create(next,packed_key);
next->lowerbound[target_state]= lowerbound;
}
add_list(next,target_state,current,source_state,added_weight,limit);

/*:93*/
#line 2748 "diag180.w"

}

/*:92*/
#line 2730 "diag180.w"

}
/*95:*/
#line 2818 "diag180.w"

table_t*swap= current;
current= next;
next= swap;
MAX_STATES= max(MAX_STATES,current->num_states);
size_t bytes= table_bytes(current)+table_bytes(next);
if(bytes> MAX_BYTES)
MAX_BYTES= bytes,MAX_RECORDS= current->num_records,MAX_KEY_BYTES= BYTES_PER_KEY;
table_finish(current);
if(!last_in_row&&(VERBOSITY>=2||(VERBOSITY==1&&wall_time()-LAST_REPORT_TIME>=PROGRESS_SECONDS)))
print_status_line("#cell",k,v,cell_number+1,num_row_cells,current,next);

/*:95*/
#line 2732 "diag180.w"


/*:91*/
#line 2493 "diag180.w"

}
/*96:*/
#line 2855 "diag180.w"

if(MERGE_MIRROR_STATES&&M_u==0&&current->num_states){
table_initialize(next,current->decided_position);
for(size_t source_state= 0;source_state<current->num_states;source_state++){
const uint8_t*packed_key= current->keys+source_state*BYTES_PER_KEY;
uint8_t mirror_key[MAX_BYTES_PER_KEY];
uint8_t mirror_crossing[MAXLAB];
coded_t coded;
expand(packed_key,current->decided_position,&coded);
int comparison= 0;
for(int a= 0,b= W-1;a<b&&comparison==0;a++,b--)
comparison= (int)coded.code[a]-MIRRORED_CODE[coded.code[b]];
if(comparison>=0){
/*98:*/
#line 2898 "diag180.w"

state_t state;
decode(&coded,&state);
int mirror_label= 0;
for(int b= W-1;b>=0;b--)
if(coded.code[b]==C_SINGLETON||coded.code[b]==C_LAST)
mirror_crossing[++mirror_label]= coded.crossing[state.label[b]];
for(int j= 1;j<=coded.num_components&&comparison==0;j++)
comparison= (int)coded.crossing[j]-(int)mirror_crossing[j];

/*:98*/
#line 2868 "diag180.w"

}
if(comparison> 0){
/*99:*/
#line 2908 "diag180.w"

coded_t mirror;
for(int a= 0;a<W;a++)
mirror.code[a]= MIRRORED_CODE[coded.code[W-1-a]];
memcpy(mirror.crossing,mirror_crossing,MAXLAB);
mirror.num_components= coded.num_components;
mirror.touches= coded.touches;
compress(&mirror,current->decided_position,mirror_key);

/*:99*/
#line 2871 "diag180.w"

packed_key= mirror_key;
}
size_t target_state= table_find(next,packed_key);
if(target_state==NO_STATE)
target_state= table_create(next,packed_key);
add_list(next,target_state,current,source_state,0,NMAX);
}
table_t*swap= current;
current= next;
next= swap;
table_finish(current);
}

/*:96*/
#line 2495 "diag180.w"

if(VERBOSITY>=1)
print_status_line("#row ",k,v,num_row_cells,num_row_cells,current,next);
/*81:*/
#line 2507 "diag180.w"

if(vspan(k-1)==W){
for(int n= 0;n<=NMAX;n++)
SQUARE[n]= add_mod(SQUARE[n],ROW_RESULTS[n]);
print_row_results(k-1);
}else if(vspan(k-1)> W){
for(int n= 0;n<=NMAX;n++)
PORTRAIT[n]= add_mod(PORTRAIT[n],ROW_RESULTS[n]);
print_row_results(k-1);
}

/*:81*/
#line 2498 "diag180.w"

}
table_free(&A);
table_free(&B);
}

/*:80*/
#line 171 "diag180.w"

/*101:*/
#line 2937 "diag180.w"

int main(int argc,char**argv)
{
/*102:*/
#line 2956 "diag180.w"

int arg= 1;
while(arg<argc&&argv[arg][0]=='-'){
if(!strcmp(argv[arg],"-m")&&arg+1<argc)
MODULUS= strtoull(argv[++arg],0,10),arg++;
else if(!strcmp(argv[arg],"-o")&&arg+1<argc)
RESULT_FILE_NAME= argv[++arg],arg++;
#line 2971 "diag180.w"
else if(!strcmp(argv[arg],"-C"))
CHECK_PACKING= 1,arg++;
else if(!strcmp(argv[arg],"-v")&&arg+1<argc)
VERBOSITY= atoi(argv[++arg]),arg++;
else if(!strcmp(argv[arg],"-p")&&arg+1<argc)
PROGRESS_SECONDS= atof(argv[++arg]),arg++;
else if(!strcmp(argv[arg],"-R"))
MERGE_MIRROR_STATES= 0,arg++;
else
break;
}
if(arg<argc&&argv[arg][0]=='-'){
fprintf(stderr,"%s: unknown option %s\n",argv[0],argv[arg]);
return 1;
}
if(arg>=argc){
fprintf(stderr,"usage: %s [options] N [Mx My umax]\n",argv[0]);
return 1;
}
NMAX= atoi(argv[arg]);
if(NMAX<1||NMAX> MAXW){
fprintf(stderr,"%s: N must be between 1 and %d\n",argv[0],MAXW);
return 1;
}
NMAX_plus_1= NMAX+1;
/*75:*/
#line 2378 "diag180.w"

PORTRAIT= calloc(NMAX_plus_1,sizeof(counter_t));
SQUARE= calloc(NMAX_plus_1,sizeof(counter_t));
ROW_RESULTS= calloc(NMAX_plus_1,sizeof(counter_t));

/*:75*/
#line 2996 "diag180.w"


/*:102*/
#line 2940 "diag180.w"

clock_t t0= clock();
START_TIME= LAST_REPORT_TIME= wall_time();
/*42:*/
#line 1288 "diag180.w"

{
uint64_t x= 0;
for(int l= 0;l<MAX_BYTES_PER_KEY;l++)
for(int b= 1;b<256;b++){
uint64_t z= (x+= 0x9E3779B97F4A7C15ull);
z= (z^(z>>30))*0xBF58476D1CE4E5B9ull;
z= (z^(z>>27))*0x94D049BB133111EBull;
HASH_BITS[l][b]= (z^(z>>31))>>32;
}
}

/*:42*/
#line 2943 "diag180.w"

/*104:*/
#line 3040 "diag180.w"

for(int stream_number= 0;stream_number<2;stream_number++){
FILE*stream= stream_number?stdout:stderr;
fprintf(stream,"# command:");
for(int word= 0;word<argc;word++)
fprintf(stream," %s",argv[word]);
fprintf(stream,"\n# diag180 built %s %s; N=%d",__DATE__,__TIME__,NMAX);
if(argc-arg==4)
fprintf(stream,", single case Mx=%s My=%s umax=%s",argv[arg+1],argv[arg+2],argv[arg+3]);
if(MODULUS)
fprintf(stream,", modulus %llu\n",(unsigned long long)MODULUS);
else
fprintf(stream,", modulus 2**"COUNTER_BITS"\n");
#line 160 "production.ch"
fprintf(stream,"# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s; pack check %s\n",
MERGE_MIRROR_STATES?"on":"off (-R)",CHECK_PACKING?"on (-C)":"off");
#line 3057 "diag180.w"
fprintf(stream,"# output: results file %s, progress -v %d -p %g s\n",
RESULT_FILE_NAME,VERBOSITY,PROGRESS_SECONDS);
fprintf(stream,"# sizes: counter %zu bytes, count record %zu bytes\n",
sizeof(counter_t),sizeof(record_t));
fflush(stream);
}

/*:104*/
#line 2944 "diag180.w"

if(argc-arg==4){
/*105:*/
#line 3066 "diag180.w"

M_x= atoi(argv[arg+1]),M_y= atoi(argv[arg+2]),umax= atoi(argv[arg+3]);
count_case();
if(MODULUS)
printf("# modulus %llu; n #(uspan<vspan) #(uspan=vspan)\n",(unsigned long long)MODULUS);
else
printf("# modulus 2**"COUNTER_BITS"; n #(uspan<vspan) #(uspan=vspan)\n");
for(int n= 1;n<=NMAX;n++)
printf("%d %s %s\n",n,counter_to_string(PORTRAIT[n]),counter_to_string(SQUARE[n]));

/*:105*/
#line 2946 "diag180.w"

}else{
/*106:*/
#line 3079 "diag180.w"


counter_t*results= calloc(3*NMAX_plus_1,sizeof(counter_t));
for(M_x= 0;M_x<=1;M_x++)
for(M_y= 0;M_y<=M_x;M_y++){


int t= M_x-M_y,center_type= M_x+M_y,twice= M_x!=M_y;
for(umax= t;2*umax-t+1<=NMAX;umax++){
count_case();
for(int n= 1;n<=NMAX;n++){
counter_t*entry= &results[center_type*NMAX_plus_1+n];
*entry= add_mod(*entry,add_mod(PORTRAIT[n],PORTRAIT[n]));
*entry= add_mod(*entry,SQUARE[n]);
if(twice){
*entry= add_mod(*entry,add_mod(PORTRAIT[n],PORTRAIT[n]));
*entry= add_mod(*entry,SQUARE[n]);
}
}
}
}

/*:106*/
#line 2948 "diag180.w"

/*107:*/
#line 3109 "diag180.w"



FILE*result_file= fopen(RESULT_FILE_NAME,"a");
if(!result_file){
perror(RESULT_FILE_NAME);
return 1;
}
char modulus_string[32];
if(MODULUS)
snprintf(modulus_string,sizeof modulus_string,"%llu",(unsigned long long)MODULUS);
else
snprintf(modulus_string,sizeof modulus_string,"2**"COUNTER_BITS);
fprintf(result_file,"# diag180 N=%d modulus %s\ntry:\n    Results\nexcept NameError:\n    Results = {}\n",NMAX,modulus_string);
static const char*type_name[3]= {"R180C","R180M","R180V"};
for(int n= 1;n<=NMAX;n++)
for(int center_type= 0;center_type<3;center_type++)
fprintf(result_file,"Results[%d, %s, '%s'] = %s\n",n,modulus_string,type_name[center_type],
counter_to_string(results[center_type*NMAX_plus_1+n]));
fclose(result_file);
printf("# modulus %s; results appended to %s\n",modulus_string,RESULT_FILE_NAME);

/*:107*/
#line 2949 "diag180.w"

/*108:*/
#line 3133 "diag180.w"

printf("# n fix(r180) R180C(center in a cell) R180M(edge midpoint) R180V(vertex)\n");
for(int n= 1;n<=NMAX;n++){
counter_t total= add_mod(add_mod(results[n],results[NMAX_plus_1+n]),results[2*NMAX_plus_1+n]);
printf("%d %s %s %s %s\n",n,counter_to_string(total),counter_to_string(results[n]),counter_to_string(results[NMAX_plus_1+n]),
counter_to_string(results[2*NMAX_plus_1+n]));
}

/*:108*/
#line 2950 "diag180.w"

}
/*109:*/
#line 3141 "diag180.w"

fprintf(stderr,"# peak memory %.1f MB (both tables; %zu count records, key %zu bytes)\n",
MAX_BYTES/1e6,MAX_RECORDS,MAX_KEY_BYTES);
fprintf(stderr,"# max states %zu, %.2fs\n",MAX_STATES,
(double)(clock()-t0)/CLOCKS_PER_SEC);

/*:109*/
#line 2952 "diag180.w"

return 0;
}

/*:101*/
#line 172 "diag180.w"


/*:2*/
