#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>



#define MAXV 500 
#define TMAX 800
#define MAXRUNTIMES 100000 


#define INPUTFILE "can__268.mtx.rnd"
#define OUTPUTFILE "ouf.txt"
#define OUT_  "OUT_t1"
#define MODE 12
#define SHAKE_MODE 1
#define SHAKE_K 0.25

struct mnode{
	int v;
	struct mnode* next;
};


struct cheak_error{
	int n;
	int VS_real;
	int VS_wrong;
};
/////////////////////////////////////////////////////


int x[MAXV] = {0};                  //解
int fx[MAXV] = {0};                 //映射表
int VS[MAXV] = {0};                 //cut成本表

int glo_x[MAXV];
int glo_fx[MAXV];
int glo_VS[MAXV];
int now_x[MAXV];
int now_fx[MAXV];
int now_VS[MAXV];



int kmax;
int c1[MAXV] = { 0 };
int c2[MAXV] = { 0 };

int vertex_num, edge_num;
int run_times = 0;
int search_times = 0;
int bug_t = 0;
int t_begin, t_end;
int imp_1=0, imp_2=0, imp_3=0;
struct mnode edge[MAXV];

struct cheak_error err;


/////////////////////////////////////////////////////


void GVNS(int* x, int *fx,int *VS,int kmax, int tmax);
void VS_updata(int* VS,int* x,int *fx, int start,int end);  //VS:原成本，x当前解
int EVS(int* x, int *VS);                                  //评价函数
void shake(int* x, int *fx, int *VS,int k);
void shake_insert(int* x, int *fx, int *VS, int k);
void constructive(int *x,int *fx,int *VS);
void input();
void NeighborhoodChange(int *x, int *fx, int *VS,
	int *x_new, int *fx_new, int *VS_new, int *k);
void insert(int *x, int *fx, int *VS,int from, int to);
int compare(int *VS1, int *VS2);
void RVNS(int *x, int *fx, int *VS, int k_max, int t_max);
void output(int*x, int *fx,int* VS);
void LocalSearch(int* x, int* fx, int* VS, int k_max);

void assignment_array(int *from, int* to, int k);
int cheak_ans(int *x, int *fx, int *VS, int mode, struct cheak_error* err);


void LocalSearch2(int* x, int* fx, int* VS, int k_max); //better than 1
void LocalSearch3(int* x, int* fx, int* VS, int k_max);
void LocalSearch4(int* x, int* fx, int* VS, int k_max);
void LocalSearch5(int* x, int* fx, int* VS, int k_max);
void LocalSearch6(int* x, int* fx, int* VS, int k_max);
void LocalSearch7(int* x, int* fx, int* VS, int k_max);
void LocalSearch8(int* x, int* fx, int* VS, int k_max);            //2领域交换改2
void LocalSearch9(int* x, int* fx, int* VS, int k_max);
void LocalSearch10(int* x, int* fx, int* VS, int k_max);
void LocalSearch11(int* x, int* fx, int* VS, int k_max);

FILE* ouf;
/////////////////////////////////////////////////////
char input_file_name[1000];
char output_file_name[1000];

int main()
{
	char tem[1000];
	unsigned int a;
	srand((unsigned)time(NULL));
	

//	printf("input filename:");
//	scanf("%s", input_file_name);

	strcpy(input_file_name, INPUTFILE);



	strcpy(output_file_name, OUT_);

	time_t timep;
	time(&timep);
	sprintf(tem,"%s", ctime(&timep));
	tem[strlen(tem)-1] = '\0';

	strcat(output_file_name, tem);
	strcat(output_file_name, input_file_name);

	

	input();
	constructive(x, fx, VS);
	kmax =(int)( SHAKE_K*vertex_num);
//	RVNS(x, fx, VS, kmax, TMAX);

	GVNS(x, fx, VS, kmax, TMAX);


	for (a = 0; a < strlen(output_file_name);a++)
	if (output_file_name[a] == ':') output_file_name[a] = '_';
	ouf = fopen(output_file_name, "w");

	cheak_ans(x, fx, VS,1,&err);

	output(x,fx,VS);



	//system("pause");
	return 0;
}

void output(int*x,int *fx ,int* VS)
{
	int a, b=0;
	
	for (a = 1; a <= vertex_num; a++)
	{
		if (VS[a] > b) b = VS[a];
	}

	fprintf(ouf, "ANS = %d\n", b);
	fprintf(ouf, "ITERATION %d TIMES\n", run_times);
	fprintf(ouf, "SEARCH %d TIMES\n", search_times);
	fprintf(ouf, "USE %f s\n", 1.0*(t_end - t_begin) / 1000);
	fprintf(ouf, "improve1 = %d ,improve2 = %d , improve3 = %d\n", imp_1, imp_2, imp_3);
	for (a = 1; a <= vertex_num; a++)
	{
		fprintf(ouf, "a=%-5dx[a]=%-8dfx[a]=%-8d  %-8d\n",a, x[a],fx[a], VS[a]);
	}
	fclose(ouf);
}

void input()
{
	int if_edge[MAXV][MAXV] = { 0 };
	char tt[1000];




	FILE *inf = fopen(input_file_name, "r");
	fgets(tt, 1000, inf);
	fscanf(inf, "%d%d%d", &vertex_num,&vertex_num, &edge_num);
	int a, b, c;
	struct mnode *p;
	for (a = 0; a <= vertex_num; a++)
	{
		edge[a].next = NULL;
		edge[a].v = a;
	}

	for (a = 0; a < edge_num; a++)
	{
		fscanf(inf, "%d%d", &b, &c);
		if (if_edge[b][c] == 0)
		{
			if_edge[b][c] = 1;
			if_edge[c][b] = 1;
			p = (struct mnode*)malloc(sizeof(struct mnode));
			p->next = edge[b].next;
			p->v = c;
			edge[b].next = p;

			p = (struct mnode*)malloc(sizeof(struct mnode));
			p->next = edge[c].next;
			p->v = b;
			edge[c].next = p;
		}
	}
	fclose(inf);
}


void VS_updata(int* VS, int* x,int *fx, int start, int end)  //VS:原成本，x当前解
{
	int incre, decre;
	int a, v_num, adj_num,max_adj;
	struct mnode *p,*q;

	if (start > end)
	{
		a = start;
		start = end;
		end = a;
	}
	for (a = start; a < end; a++)
	{
		incre = 0;   
		decre = 0;
		v_num = x[a];
		p = edge[v_num].next;

		while (p)
		{
			adj_num=fx[p->v];
			if (adj_num>a)
			{
				incre = 1;
			}
			if (adj_num < a)
			{
				q = edge[x[adj_num]].next;
				max_adj = 0;
				while (q)
				{
					if (fx[q->v] > max_adj) max_adj = fx[q->v];
				   	q = q->next;
				}
				if (max_adj == a)
					decre--;
			}
	    	p = p->next;
		}
		if (a == 1) VS[a] = 1;
		else
		{
			VS[a] = VS[a - 1] + incre + decre;
		}
		
	}
	return;
}

int compare(int *VS1, int *VS2)
{
	int a, maxc;
	
	maxc = 0;
	for (a = 0; a < vertex_num; a++)
	{
		c1[a] = 0;
		c2[a] = 0;
	}
	for (a = 0; a < vertex_num; a++)
	{
		c1[VS1[a]]++;
 		c2[VS2[a]]++;
		if (VS1[a] > maxc) maxc = VS1[a];
		if (VS2[a] > maxc) maxc = VS2[a];
	}
	for (a = maxc; a >0; a--)
	{
		if (c1[a] > c2[a]) return 2;          //VS2 better
		if (c1[a] < c2[a]) return 1;          //VS1 better
	}
	return 3;              //the same
}


void shake(int* x, int *fx, int *VS, int k)
{
	int a, b, c, d;
	for (a = 0; a < k; a++)
	{
		b = (rand() % vertex_num) + 1 ;
		c = (rand() % vertex_num) + 1;
		fx[x[b]] = c;
		fx[x[c]] = b;
		d = x[b];
		x[b] = x[c];
		x[c] = d;

		VS_updata(VS, x, fx, b, c);
	}
	return;
}

void assignment_array(int *from, int* to,int k)
{
	int a;
	for (a = 0; a <= k; a++) to[a] = from[a];
	return;
}

void NeighborhoodChange(int *x, int *fx,int *VS,
	int *x_new, int *fx_new,int *VS_new, int *k)
{
	if (compare(VS, VS_new) == 2)
	{
		assignment_array(VS_new, VS, vertex_num);
		assignment_array(x_new,x,vertex_num);
		assignment_array(fx_new,fx,vertex_num);
		*k = 1;
	}
	else
	{
		(*k)++;
	}
	return;
}

void shake_insert(int* x, int *fx, int *VS, int k)
{
	int a, b, c;
	int from, to;
	for (a = 0; a < k; a++)
	{
		b = (rand() % vertex_num) + 1;
		c = (rand() % vertex_num) + 1;
		from = b;
		to = c;
	//	insert(x, fx, VS, b, c);
		{
			int a, b;
			if (from < to)
			{
				a = x[from];
				for (b = from; b <= to - 1; b++)
				{
					x[b] = x[b + 1];
					fx[x[b]] = b;
				}
				x[to] = a;
				fx[x[to]] = to;
			}
			else
			{
				a = x[from];
				for (b = from; b >= to + 1; b--)
				{
					x[b] = x[b - 1];
					fx[x[b]] = b;
				}
				x[to] = a;
				fx[x[to]] = to;
			}
		}
	
	}
	VS_updata(VS, x, fx, 1, vertex_num);
	return;
}

void insert(int *x, int *fx, int *VS, int from, int to)
{
	int a, b;

	if (from < to)
	{
		a = x[from];
		for (b = from; b<=to-1 ; b++)
		{
			x[b] = x[b + 1];
			fx[x[b]] = b;
		}
		x[to] = a;
		fx[x[to]] = to;
	}
	else
	{
		a = x[from];
		for (b = from; b >= to+1; b--)
		{
			x[b] = x[b - 1];
			fx[x[b]] = b;
		}
		x[to] = a;
		fx[x[to]] = to;
	}
	
	VS_updata(VS, x, fx, from, to);
	return;
}

void swap(int *x, int *fx, int *VS, int from, int to)
{
	int a;
	a = x[from];
	x[from] = x[to];
	x[to] = a;
	fx[x[from]] = from;
	fx[x[to]] = to;

	VS_updata(VS, x, fx, from, to);
	return;
}



void constructive(int *x, int *fx, int *VS)
{
	int BFS_que[MAXV] = { 0 };
	int if_search[MAXV] = { 0 };
	int ans[MAXV] = { 0 };
	int depth_tem[MAXV] = { 0 };
	int tem_x[MAXV] = { 0 }, tem_VS[MAXV] = { 0 }, tem_fx[MAXV] = { 0 };
	int max_x;
	int last, now;
	int a, b, c;
	int depth=0,tem_max_depth;
	struct mnode *p;

	for (a = 1; a <= vertex_num; a++)
	{
		now = 1;
		last = 1;
		BFS_que[1] = a;
		tem_max_depth = 1;
		for (b = 1; b <= vertex_num; b++)
		{
			if_search[b] = 0;
			depth_tem[b] = 0;
		}
		if_search[a] = 1;
		depth_tem[a] = 1;

		while (now <= last)
		{
			p = edge[BFS_que[now]].next;
			while (p)
			{
				if (if_search[p->v] == 0)
				{
					++last;
					BFS_que[last] = p->v;
					depth_tem[p->v] = depth_tem[BFS_que[now]]+1;
					if (depth_tem[p->v] > tem_max_depth)
						tem_max_depth = depth_tem[p->v];
					if_search[p->v]=1;
				}
				p = p->next;
			}
			now++;
		}
		if (tem_max_depth > depth)
		{
			depth = tem_max_depth;
			assignment_array(BFS_que, ans, vertex_num);
		}
	}

	x[1] = ans[1];
	fx[ans[1]] = 1;
	VS[1] = 0;

	

	for (a = 2; a <= vertex_num; a++)
	{
		max_x = 0;
		for (b = a; b >= 2; b--)
		{
			x[b] = x[b - 1];
			VS[b] = VS[b - 1];
			fx[x[b]] = b;
		}
		x[1] = ans[a];
		VS[1] = 0;
		fx[ans[a]] = 1;

		p = edge[ans[a]].next;
		while (p)
		{
			//////fx[k]=0 当k未被初始化
			if (fx[p->v] > max_x)
				max_x = fx[p->v];
			p = p->next;
		}
		for (b = 1; b < max_x; b++)
			VS[b]++;

		///插入到开头，更新VS，调换位置寻找最优VS
		assignment_array(x, tem_x, a);
		assignment_array(VS, tem_VS, a);
		assignment_array(fx, tem_fx, vertex_num);
		
		///两两交换判定是否更优

		for (b = 1; b < a; b++)
		{
			c = tem_x[b];
			tem_x[b] = tem_x[b + 1];
			tem_x[b + 1] = c;
			tem_fx[tem_x[b]] = b;
			tem_fx[tem_x[b + 1]] = b + 1;
			VS_updata(tem_VS, tem_x, tem_fx, b, b+1);
			c = compare(tem_VS, VS);
			if (c == 1)
			{
				assignment_array(tem_x, x, a);
				assignment_array(tem_VS, VS, a);
				assignment_array(tem_fx, fx, vertex_num);
			}

	//	if (cheak_ans(tem_x, tem_fx, tem_VS, 0, &err) == -1)
		//	{
	      //	printf("BUGBUGBUGBUGBUGBUGBUGBUG\n");
		//	}


		}
		
		///生成初始解
	}

}




void RVNS(int *x, int *fx, int *VS, int k_max, int t_max)
{
	int k;
	int tem_x[MAXV] = { 0 }, tem_VS[MAXV] = { 0 }, tem_fx[MAXV] = {0};
	t_begin = clock();
	do
	{
		k = 1;
		do
		{
			bug_t++;
			assignment_array(x, tem_x, vertex_num);
			assignment_array(fx, tem_fx, vertex_num);
			assignment_array(VS, tem_VS, vertex_num);

		//	if (cheak_ans(tem_x, tem_fx, tem_VS, 0, &err) == -1)
	//		{
//				printf("BUGBUGBUGBUGBUGBUGBUGBUG\n");
	//		}

			shake(tem_x, tem_fx, tem_VS, k);

//			if (cheak_ans(tem_x, tem_fx, tem_VS, 0, &err) == -1)
//			{
//				printf("BUGBUGBUGBUGBUGBUGBUGBUG\n");
	//		}

			NeighborhoodChange(x, fx, VS, tem_x, tem_fx, tem_VS, &k);


		} while (k <= k_max);
		t_end = clock();
	}while (1.0*(t_end - t_begin) / 1000 < t_max);
	return;
}


int cheak_ans(int *x,int *fx,int *VS,int mode,struct cheak_error* err)
{
	int a, b, c;
	int cheak[MAXV] = {0};
	struct mnode *p;
	for (a = 1; a <= vertex_num; a++)
	{
		c = 0;
		p = edge[x[a]].next;
		while (p) 
		{ 
			if (fx[p->v] > c) c = fx[p->v]; 
			p = p->next; 
		} 
		for (b = a; b < c; b++) cheak[b]++; 
	} 
	for (a = 1; a <= vertex_num; a++) 
	{ 
		if (cheak[a] != VS[a]) 
		{   
			if (mode==1) 
		    	fprintf(ouf,"ANS CHEAK WRONG in %d  %d!=%d\n", a, cheak[a], VS[a]); 
			err->n = a; 
			err->VS_real = cheak[a]; 
			err->VS_wrong = VS[a]; 
			return -1;  
		} 
	} 
	if (mode==1) 
	    fprintf(ouf,"ANS CHEAK SUCCESS\n"); 
	return 1; 
} 

int loc_x[MAXV];
int loc_fx[MAXV];
int loc_VS[MAXV];


void LocalSearch(int* x, int* fx, int* VS, int k_max)
{
	int improvement = 1;
	int a, b, c ;
	while (improvement)
	{
		improvement = 0;
		for (a = 1; a <= vertex_num; a++)
		{
			for (b = 1; b <= vertex_num; b++)
			{
				if (a != b)
				{
					++search_times;
					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
				//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}

				}
			}
		}

	}
	return;
}

void LocalSearch2(int* x, int* fx, int* VS, int k_max)
{
	int improvement = 1;
	int a, b, c, d;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	while (improvement)
	{
		improvement = 0;

		for (a = 1; a <= vertex_num; a++)
		{
			for (b = 1; b <= vertex_num; b++)
			{
				if (a != b)
				{
					++search_times;
				//	assignment_array(x, loc_x, vertex_num);
				//	assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					} 
					else
					{
						d = loc_x[a];
						loc_x[a] = loc_x[b];
						loc_x[b] = d;
						loc_fx[loc_x[a]] = a;
						loc_fx[loc_x[b]] = b;
					}

				}
			}
		}
	}
	return;
}

void LocalSearch3(int* x, int* fx, int* VS, int k_max)
{
	int improvement = 1;
	int a, b, c;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	while (improvement)
	{
		improvement = 0;

		for (a = 1; a <= vertex_num; a++)
		{
			m1 = 9999;  m2 = 9999;
			p = edge[x[a]].next;
			while (p)
			{
				if (fx[p->v] < m1)
				{
					m2 = m1;
					m1 = fx[p->v];
				}
				else
				{
					if (fx[p->v] < m2)
					{
						m2 = fx[p->v];
					}
				}
				p = p->next;
			}
			
			if (m1 == vertex_num) b = vertex_num;
			else
				b = m1 + 1 + rand() % (vertex_num - m1);

			if (b > vertex_num)b = vertex_num;
			if (a != b)
			{
				++search_times;

				assignment_array(x, loc_x, vertex_num);
				assignment_array(fx, loc_fx, vertex_num);
				assignment_array(VS, loc_VS, vertex_num);
			//	insert(loc_x, loc_fx, loc_VS, a, b);
				swap(loc_x, loc_fx, loc_VS, a, b);
				c = compare(VS, loc_VS);
				if (c == 2)
				{
					improvement = 1;
					assignment_array(loc_x, x, vertex_num);
					assignment_array(loc_fx, fx, vertex_num);
					assignment_array(loc_VS, VS, vertex_num);
				}
			}
			
		}
	}
	return;
}


void LocalSearch4(int* x, int* fx, int* VS, int k_max)
{
	int improvement = 1;
	int a, b, c;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{
		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
				}

			}
		}
		if (!improvement)
		{
			for (a = 1; a <= vertex_num; a++)
			{
				for (b = 1; b <= vertex_num; b++)
				{
					if (a != b)
					{
						++search_times;
						assignment_array(x, loc_x, vertex_num);
						assignment_array(fx, loc_fx, vertex_num);
						assignment_array(VS, loc_VS, vertex_num);
						//	insert(loc_x, loc_fx, loc_VS, a, b);
						swap(loc_x, loc_fx, loc_VS, a, b);
						c = compare(VS, loc_VS);
						if (c == 2)
						{
							improvement = 1;
							assignment_array(loc_x, x, vertex_num);
							assignment_array(loc_fx, fx, vertex_num);
							assignment_array(loc_VS, VS, vertex_num);
						}

					}
				}
			}
		}
	} while (improvement);
	return;
}


void LocalSearch5(int* x, int* fx, int* VS, int k_max)            //2领域交换
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
				}

			}
		}
		if (!improvement)
		{

			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
						//	insert(loc_x, loc_fx, loc_VS, a, b);
					 		swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
							}

						}
					}
				}
			}
		}
		
	} while (improvement);
	return;
}


void LocalSearch6(int* x, int* fx, int* VS, int k_max)         //3领域
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						++imp_1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
				}

			}
		}
		if (!improvement)
		{
			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							//	insert(loc_x, loc_fx, loc_VS, a, b);
							swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								++imp_2;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
							}

						}
					}
				}
				if (!improvement)
				{
					for (a = 1; a <= vertex_num; a++)
					{
						for (b = 1; b <= vertex_num; b++)
						{
							if (a != b)
							{
								++search_times;
								assignment_array(x, loc_x, vertex_num);
								assignment_array(fx, loc_fx, vertex_num);
								assignment_array(VS, loc_VS, vertex_num);
								insert(loc_x, loc_fx, loc_VS, a, b);
							//  swap(loc_x, loc_fx, loc_VS, a, b);
								c = compare(VS, loc_VS);
								if (c == 2)
								{
									improvement = 1;
									++imp_3;
									assignment_array(loc_x, x, vertex_num);
									assignment_array(loc_fx, fx, vertex_num);
									assignment_array(loc_VS, VS, vertex_num);
								}

							}
						}
					}
				}
			}
		}
	

	} while (improvement);
	return;
}

void LocalSearch7(int* x, int* fx, int* VS, int k_max)            //2领域交换改1
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	int ls_x[MAXV], ls_fx[MAXV], ls_VS[MAXV];
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);

	

	do
	{
		

		while (improvement)
		{
			improvement = 0;
			assignment_array(x, ls_x, vertex_num);
			assignment_array(fx, ls_fx, vertex_num);
			assignment_array(VS, ls_VS, vertex_num);


			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);

					c = compare(ls_VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, ls_x, vertex_num);
						assignment_array(loc_fx, ls_fx, vertex_num);
						assignment_array(loc_VS, ls_VS, vertex_num);
					}
				}

			}

			
			assignment_array(ls_x, x, vertex_num);
			assignment_array(ls_fx, fx, vertex_num);
			assignment_array(ls_VS, VS, vertex_num);
		
		}
		if (!improvement)
		{

			d = compare(ls_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							//	insert(loc_x, loc_fx, loc_VS, a, b);
							swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(ls_VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								assignment_array(loc_x, ls_x, vertex_num);
								assignment_array(loc_fx, ls_fx, vertex_num);
								assignment_array(loc_VS, ls_VS, vertex_num);
							}

						}
					}
				}
				assignment_array(ls_x, x, vertex_num);
				assignment_array(ls_fx, fx, vertex_num);
				assignment_array(ls_VS, VS, vertex_num);
			}
		}


	} while (improvement);
	return;
}

void LocalSearch8(int* x, int* fx, int* VS, int k_max)            //2临域交换改2
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	int ifbreak;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
						break;
					}
				}

			}
		}
		if (!improvement)
		{
			ifbreak = 0;
			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							//	insert(loc_x, loc_fx, loc_VS, a, b);
							swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
								ifbreak = 1;
								break;
							}

						}
					}
					if (ifbreak) break;
				}
			}
		}

	} while (improvement);
	return;
}


void LocalSearch9(int* x, int* fx, int* VS, int k_max)         //3领域,n临域采用insert方法
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

			

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;

					assignment_array(x, loc_x, vertex_num);
					assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					insert(loc_x, loc_fx, loc_VS, a, b);
					//swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
						++imp_1;
					}
				}

			}
		}
		if (!improvement)
		{
			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							//	insert(loc_x, loc_fx, loc_VS, a, b);
							swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
								++imp_2;
							}

						}
					}
				}
				if (!improvement)
				{
					for (a = 1; a <= vertex_num; a++)
					{
						for (b = 1; b <= vertex_num; b++)
						{
							if (a != b)
							{
								++search_times;
								assignment_array(x, loc_x, vertex_num);
								assignment_array(fx, loc_fx, vertex_num);
								assignment_array(VS, loc_VS, vertex_num);
								insert(loc_x, loc_fx, loc_VS, a, b);
								//  swap(loc_x, loc_fx, loc_VS, a, b);
								c = compare(VS, loc_VS);
								if (c == 2)
								{
									improvement = 1;
									assignment_array(loc_x, x, vertex_num);
									assignment_array(loc_fx, fx, vertex_num);
									assignment_array(loc_VS, VS, vertex_num);
									++imp_3;
								}

							}
						}
					}
				}
			}
		}


	} while (improvement);
	return;
}


void LocalSearch10(int* x, int* fx, int* VS, int k_max)         //3领域，n^2临域先插入后交换，优化实现
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;
					//	assignment_array(x, loc_x, vertex_num);
					//	assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;

						++imp_1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
					else
					{
						d = loc_x[a];
						loc_x[a] = loc_x[b];
						loc_x[b] = d;
						loc_fx[loc_x[a]] = a;
						loc_fx[loc_x[b]] = b;
					}

				}

			}
		}
		if (!improvement)
		{
			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							insert(loc_x, loc_fx, loc_VS, a, b);
							//swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								++imp_2;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
							}

						}
					}
				}
				if (!improvement)
				{
					for (a = 1; a <= vertex_num; a++)
					{
						for (b = 1; b <= vertex_num; b++)
						{
							if (a != b)
							{
								++search_times;
								assignment_array(x, loc_x, vertex_num);
								assignment_array(fx, loc_fx, vertex_num);
								assignment_array(VS, loc_VS, vertex_num);
								//insert(loc_x, loc_fx, loc_VS, a, b);
								swap(loc_x, loc_fx, loc_VS, a, b);
								c = compare(VS, loc_VS);
								if (c == 2)
								{
									improvement = 1;
									++imp_3;
									assignment_array(loc_x, x, vertex_num);
									assignment_array(loc_fx, fx, vertex_num);
									assignment_array(loc_VS, VS, vertex_num);
								}

							}
						}
					}
				}
			}
		}


	} while (improvement);
	return;
}

void LocalSearch11(int* x, int* fx, int* VS, int k_max)         //3领域，优化赋值
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d, e;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;
					//	assignment_array(x, loc_x, vertex_num);
					//	assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;

						++imp_1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
					else
					{
						d = loc_x[a];
						loc_x[a] = loc_x[b];
						loc_x[b] = d;
						loc_fx[loc_x[a]] = a;
						loc_fx[loc_x[b]] = b;
					}

				}

			}
		}
		if (!improvement)
		{
			d = compare(loc_VS, glo_VS);
			if (d != 2)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
						//	assignment_array(x, loc_x, vertex_num);  
						//	assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
						//	insert(loc_x, loc_fx, loc_VS, a, b); 
							swap(loc_x, loc_fx, loc_VS, a, b);
							c = compare(VS, loc_VS); 
							if (c == 2) 
							{ 
								improvement = 1;   
								++imp_2;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num); 
								assignment_array(loc_VS, VS, vertex_num); 
							}
							else 
							{
								e = loc_x[a]; 
								loc_x[a] = loc_x[b];
								loc_x[b] = e;   
								loc_fx[loc_x[a]] = a;
								loc_fx[loc_x[b]] = b;
							}  
						}
					} 
				}
				if (!improvement)  
				{
					for (a = 1; a <= vertex_num; a++) 
					{ 
						for (b = 1; b <= vertex_num; b++)
						{ 
							if (a != b) 
							{ 
								++search_times;
								assignment_array(x, loc_x, vertex_num); 
								assignment_array(fx, loc_fx, vertex_num); 
								assignment_array(VS, loc_VS, vertex_num);
								insert(loc_x, loc_fx, loc_VS, a, b);
							//	swap(loc_x, loc_fx, loc_VS, a, b);  
								c = compare(VS, loc_VS); 
								if (c == 2)
								{ 
									improvement = 1;   
									++imp_3;
									assignment_array(loc_x, x, vertex_num);
									assignment_array(loc_fx, fx, vertex_num);
									assignment_array(loc_VS, VS, vertex_num);
								}

							}
						}
					}
				}
			}
		} 


	} while (improvement);
	return;
}


void LocalSearch12(int* x, int* fx, int* VS, int k_max)         //3领域，优化赋值 ,领域变换不考虑全局解
{
	int improvement = 1;
	int improvement_glo = 0;
	int a, b, c, d, e;
	int m1, m2;
	struct mnode *p;
	assignment_array(x, loc_x, vertex_num);
	assignment_array(fx, loc_fx, vertex_num);
	assignment_array(VS, loc_VS, vertex_num);
	do
	{

		while (improvement)
		{
			improvement = 0;

			for (a = 1; a <= vertex_num; a++)
			{
				m1 = 9999;  m2 = 9999;
				p = edge[x[a]].next;
				while (p)
				{
					if (fx[p->v] < m1)
					{
						m2 = m1;
						m1 = fx[p->v];
					}
					else
					{
						if (fx[p->v] < m2)
						{
							m2 = fx[p->v];
						}
					}
					p = p->next;
				}

				if (m1 == vertex_num) b = vertex_num;
				else
					b = m1 + 1 + rand() % (vertex_num - m1);

				if (b > vertex_num)b = vertex_num;
				if (a != b)
				{
					++search_times;
					//	assignment_array(x, loc_x, vertex_num);
					//	assignment_array(fx, loc_fx, vertex_num);
					assignment_array(VS, loc_VS, vertex_num);
					//	insert(loc_x, loc_fx, loc_VS, a, b);
					swap(loc_x, loc_fx, loc_VS, a, b);
					c = compare(VS, loc_VS);
					if (c == 2)
					{
						improvement = 1;

						++imp_1;
						assignment_array(loc_x, x, vertex_num);
						assignment_array(loc_fx, fx, vertex_num);
						assignment_array(loc_VS, VS, vertex_num);
					}
					else
					{
						d = loc_x[a];
						loc_x[a] = loc_x[b];
						loc_x[b] = d;
						loc_fx[loc_x[a]] = a;
						loc_fx[loc_x[b]] = b;
					}

				}

			}
		}
		if (!improvement)
		{
			for (a = 1; a <= vertex_num; a++)
			{
				for (b = 1; b <= vertex_num; b++)
				{
					if (a != b)
					{
						++search_times;
						//	assignment_array(x, loc_x, vertex_num);  
						//	assignment_array(fx, loc_fx, vertex_num);
						assignment_array(VS, loc_VS, vertex_num);
						//	insert(loc_x, loc_fx, loc_VS, a, b); 
						swap(loc_x, loc_fx, loc_VS, a, b);
						c = compare(VS, loc_VS);
						if (c == 2)
						{
							improvement = 1;
							++imp_2;
							assignment_array(loc_x, x, vertex_num);
							assignment_array(loc_fx, fx, vertex_num);
							assignment_array(loc_VS, VS, vertex_num);
						}
						else
						{
							e = loc_x[a];
							loc_x[a] = loc_x[b];
							loc_x[b] = e;
							loc_fx[loc_x[a]] = a;
							loc_fx[loc_x[b]] = b;
						}
					}
				}
			}
			if (!improvement)
			{
				for (a = 1; a <= vertex_num; a++)
				{
					for (b = 1; b <= vertex_num; b++)
					{
						if (a != b)
						{
							++search_times;
							assignment_array(x, loc_x, vertex_num);
							assignment_array(fx, loc_fx, vertex_num);
							assignment_array(VS, loc_VS, vertex_num);
							insert(loc_x, loc_fx, loc_VS, a, b);
							//	swap(loc_x, loc_fx, loc_VS, a, b);  
							c = compare(VS, loc_VS);
							if (c == 2)
							{
								improvement = 1;
								++imp_3;
								assignment_array(loc_x, x, vertex_num);
								assignment_array(loc_fx, fx, vertex_num);
								assignment_array(loc_VS, VS, vertex_num);
							}

						}
					}
				}
			}
			
		}


	} while (improvement);
	return;
}



void GVNS(int* x, int *fx, int *VS, int shakek, int tmax)
{
	int a;


	assignment_array(x, glo_x,vertex_num);
	assignment_array(fx, glo_fx, vertex_num);
	assignment_array(VS, glo_VS, vertex_num);

	assignment_array(x, now_x, vertex_num);
	assignment_array(fx, now_fx, vertex_num);
	assignment_array(VS, now_VS, vertex_num);

	t_begin = clock(); 
	do
	{
		assignment_array(glo_x, now_x, vertex_num);
		assignment_array(glo_fx, now_fx, vertex_num);
		assignment_array(glo_VS, now_VS, vertex_num);

		if (SHAKE_MODE == 1)
		{
			shake(now_x, now_fx, now_VS, shakek);
		}
		else
		{
			shake_insert(now_x, now_fx, now_VS, shakek);
		}

		++run_times;

		switch (MODE)
		{
		case 1:
		{
				  LocalSearch(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 2:
		{
				  LocalSearch2(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 3:
		{
				  LocalSearch3(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 4:
		{
				  LocalSearch4(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 5:
		{
				  LocalSearch5(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 6:
		{
				  LocalSearch6(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 7:
		{
				  LocalSearch7(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 8:
		{
				  LocalSearch8(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 9:
		{
				  LocalSearch9(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 10:
		{
				  LocalSearch10(now_x, now_fx, now_VS, 0);
				  break;
		}
		case 11:
		{
				   LocalSearch11(now_x, now_fx, now_VS, 0);
				   break;
		}
		case 12:
		{
				   LocalSearch12(now_x, now_fx, now_VS, 0);
				   break;
		}
		
		default:
			break;
		}
		//LocalSearch(now_x, now_fx, now_VS, 0);
	//	LocalSearch2(now_x, now_fx, now_VS, 0);
	//	LocalSearch3(now_x, now_fx, now_VS, 0);
	//	LocalSearch4(now_x, now_fx, now_VS, 0);
	//	LocalSearch5(now_x, now_fx, now_VS, 0);
	//	LocalSearch6(now_x, now_fx, now_VS, 0);
	//	LocalSearch7(now_x, now_fx, now_VS, 0);
	//	LocalSearch8(now_x, now_fx, now_VS, 0);
		a = compare(glo_VS, now_VS);
		if (a == 2)
		{
			assignment_array(now_x, glo_x, vertex_num);
			assignment_array(now_fx, glo_fx, vertex_num);
			assignment_array(now_VS, glo_VS, vertex_num);
		}
		t_end = clock();
	} while ((1.0*(t_end - t_begin) / 1000 < tmax) && (run_times < MAXRUNTIMES));


	assignment_array(glo_x, x, vertex_num);
	assignment_array(glo_fx, fx, vertex_num);
	assignment_array(glo_VS, VS, vertex_num);

	return;
}