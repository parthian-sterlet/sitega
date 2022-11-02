#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
//#include  <conio.h>
#include  <math.h>
#include  <time.h>
#include <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define MOTLEN 12 //max LPD length
#define MEGE 40//population size 1st stage
#define ELIT 20//population size 2nd stage
#define NMUT 3
#define NREC 5
#define POPSIZE 110
#define CENT 100
#define GA_EXIT 0.05

double uw[POPSIZE][POPSIZE];

struct ss {
	int num;
	char oli[3];
}s[16];

int compare_qq(const void* X1, const void* X2)//increase
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return 1;
	if (X < 0)return -1;
	return 0;
}
int compare_qq2(const void* X1, const void* X2)//decrease
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return -1;
	if (X < 0)return 1;
	return 0;
}
struct qbs {
	double q;//best score
	int n;// 1 forgr 0 backgr
};
int compare_qbs(const void* X1, const void* X2)//decrease
{
	struct qbs* S1 = (struct qbs*)X1;
	struct qbs* S2 = (struct qbs*)X2;
	if (S1->q - S2->q > 0)return -1;
	if (S1->q - S2->q < 0)return 1;
	if (S1->n - S2->n > 0)return 1;
	if (S1->n - S2->n < 0)return -1;
	return 0;
}
struct uno {
	int sta;
	int end;
	int num;
	void get_copy(uno* a);
	void print_all(void);
};
void uno::get_copy(uno* a)
{
	a->num = num;
	a->sta = sta;
	a->end = end;
};
void uno::print_all(void)
{
	printf("[%d;%d]%s\t", sta, end, s[num].oli);
}
struct town {
	uno tot[POPSIZE];
	int deg[16];
	int size;
	int* pos;// pozicii na4al okon
	//int *posc;// pozicii na4al okon
	int* ori;// DNA strand 0,1
	double fit;
	double mah;
	double fpr;
	int sm;// success mut
	int sr;// success rec
	int tm;// try mut
	int tr;// try rec
	int odg[MOTLEN + 1];
	void get_copy(town* a, int nseq, int reg_max);
	//void init_rand(int nseq, int *len, int oln, int rsize, int reg_max);
	void init_rand_hoxa(int nseq, int oln, int rsize, int reg_max, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	//void init_rand_part(int nseq, int *len, int oln, int *xporti, int nind);
	void init_rand_part_hoxa(int nseq, int* xporti, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	int init_add(uno last);
	int init_add_split(void);
	void init_zero(int olen);
	void Mix(int* a, int* b);
	void BigMix(int* d, int ln);
	void swap(int n1, int n2);
	int order(int n);
	void print_all(int reg_max, int nseq);
	void print_sta(int reg_max);
	void sort_all(void);
	int sum(int j);
	void fprint_all(char* file, char* add);
	void fprint_allfi(char* file, char* add, int len, double sd, double c0, double* buf);
	void fprint_allfi_mat(char* file, char* add, char* name, int len, double c0, double* buf, int iter, int size_start, int olen_min);
	int check(int min, int max);
	int mem_in(int nseq);
	void mem_out(void);
}**pop, det1, det2[2];
int compare_tot(const void* X1, const void* X2)
{
	struct uno* S1 = (struct uno*)X1;
	struct uno* S2 = (struct uno*)X2;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	if (S1->sta - S2->sta > 0)return 1;
	if (S1->sta - S2->sta < 0)return -1;
	return 0;
}
int compare_pop(const void* X1, const void* X2)
{
	struct town* S1 = (struct town*)X1;
	struct town* S2 = (struct town*)X2;
	if (S1->fit - S2->fit > 0)return -1;
	if (S1->fit - S2->fit < 0)return 1;
	return 0;
}
void town::init_zero(int reg_max)
{
	int k;
	for (k = 0; k < 16; k++)deg[k] = 0;
	for (k = 0; k < reg_max; k++)odg[k] = 0;
	odg[reg_max] = -1;
	size = 0;
	fit = mah = 0;
	fpr = 1;
}
int town::init_add(uno last)
{
	int k;
	int dtype = last.num;
	int sta = last.sta;
	int end = last.end;
	for (k = 0; k < size; k++)
	{
		if (tot[k].num == dtype)
		{
			int sta1 = tot[k].sta;
			int end1 = tot[k].end;
			if (sta >= sta1 && sta <= end1)return -1;
			if (end >= sta1 && end <= end1)return -1;
			if (sta1 >= sta && sta1 <= end)return -1;
			if (end1 >= sta && end1 <= end)return -1;
		}
	}
	tot[size].sta = sta;
	tot[size].end = end;
	tot[size].num = dtype;
	deg[dtype]++;
	odg[end - sta]++;
	size++;
	return 1;
}
int town::init_add_split(void)
{
	int k, r, r1, x, end1;
	do
	{
		r = rand() % size;
	} while (tot[r].end - tot[r].sta == 0);
	for (k = size - 1; k > r; k--)
	{
		tot[k].get_copy(&tot[k + 1]);
	}
	x = tot[r].end - tot[r].sta;
	end1 = rand() % x;
	r1 = r + 1;
	tot[r1].sta = tot[r].sta + end1 + 1;
	tot[r1].end = tot[r].end;
	tot[r1].num = tot[r].num;
	tot[r].end = tot[r].sta + end1;
	deg[tot[r].num]++;
	odg[x]--;
	odg[tot[r1].end - tot[r1].sta]++;
	odg[tot[r].end - tot[r].sta]++;
	size++;
	return 1;
}
void town::sort_all(void)
{
	qsort((void*)tot, size, sizeof(tot[0]), compare_tot);
}
int town::check(int min, int max)
{
	int i, j;
	int ods[MOTLEN];
	int max1 = max - 1;
	for (i = 0; i < max1; i++)ods[i] = 0;
	for (i = 0; i < size; i++)
	{
		ods[tot[i].end - tot[i].sta]++;
	}
	for (i = 0; i < max1; i++)
	{
		if (ods[i] != odg[i])
		{
			printf("odg(%d) %d actual %d", i, odg[i], ods[i]);
			return -1;
		}
	}
	int sum = 0;
	for (i = 0; i < max; i++)
	{
		if (odg[i] < 0) { printf("odg %d = %d", i, odg[i]); return -1; }
		sum += odg[i];
	}
	if (sum != size)
	{
		printf("wrong sum odg = %d, size = %d", sum, size);
		return -1;
	}
	for (i = 0; i < size; i++)
	{
		int len = tot[i].end - tot[i].sta;
		if (len < min || len >= max)
		{
			tot[i].print_all();
			printf("\n");
			return -1;
		}
	}
	for (i = 0; i < size; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			if (tot[i].num == tot[j].num)
			{
				if (tot[i].end <= tot[j].sta)continue;
				if (tot[j].end <= tot[i].sta)continue;
				int cat = 0;
				if (tot[i].sta >= tot[j].sta && tot[j].sta <= tot[i].end)cat = 1;
				if (tot[j].sta >= tot[i].sta && tot[i].sta <= tot[j].end)cat = 1;
				if (tot[i].sta >= tot[j].end && tot[j].end <= tot[i].end)cat = 1;
				if (tot[j].sta >= tot[i].end && tot[i].end <= tot[j].end)cat = 1;
				if (cat == 1)
				{
					printf("1st %d\t", i + 1);
					tot[j].print_all();
					printf("\n");
					printf("2nd %d\t", j + 1);
					tot[i].print_all();
					printf("\n");
					return -1;
				}
			}
		}
	}
	return 1;
}
void town::fprint_allfi_mat(char* file, char* add, char* name, int len, double c0, double* buf, int iter, int size_start, int olen_min)
{
	int i;
	FILE* out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);

	if (iter == 0 && (olen_min == len && size_start == size))
	{
		if ((out = fopen(file_out, "wt")) == NULL)
		{
			printf("Ouput file can't be opened!\n");
			exit(1);
		}
	}
	else
	{
		if ((out = fopen(file_out, "at")) == NULL)
		{
			printf("Ouput file can't be opened!\n");
			exit(1);
		}
	}
	fprintf(out, "%s_%d\n", name, iter);
	fprintf(out, "%d\tLPD count\n", size);
	fprintf(out, "%d\tModel length\n", len);
	fprintf(out, "%.12f\tCoefficient\n", c0);
	for (i = 0; i < size; i++)
	{
		fprintf(out, "%d\t%d\t%.12f\t", tot[i].sta, tot[i].end, buf[i]);
		fprintf(out, "%d\t%s\n", tot[i].num, s[tot[i].num].oli);
	}
	fclose(out);
}
struct town_ext {
	double c0;
	double buf[POPSIZE];
	void get_copy(double c1, double* b1, int size);
} pop_ext;
void town_ext::get_copy(double c1, double* b1, int size)
{
	c0 = c1;
	int i;
	for (i = 0; i < size; i++)buf[i] = b1[i];
}
void MixI(int* a, int* b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixI(int* d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r < len - 1; r++)
	{
		MixI(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void MixC(char* a, char* b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixC(char* d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r < len - 1; r++)
	{
		MixC(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void town::print_all(int reg_max, int nseq)
{
	int i;
	char strand[] = "+-";
	printf("M %f FP %f F %f S %d\t", mah, fpr, fit, size);
	for (i = 0; i < 16; i++)printf("%d ", deg[i]); printf("\tLEN ");
	for (i = 0; i < reg_max; i++)printf("%d ", odg[i]); printf("\t");
	int size1 = size - 1;
	printf(" ");
	for (i = 0; i < size; i++)
	{
		tot[i].print_all();
		if (i == size1)printf("\n");
		// else if(tot[i + 1].num != tot[i].num)printf("\n");
		//	 printf(" ");
	}
	// printf("\n");

	for (i = 0; i < nseq; i++)
	{
		printf("%2d", pos[i]);
		printf("%c", strand[ori[i]]);
		printf(" ");
	}
	printf("\n");
}
void town::print_sta(int reg_max)
{
	int i;
	printf("M %f H %f F %f\t", mah, fpr, fit);
	for (i = 0; i < 16; i++)printf("%2d", deg[i]); printf("\t");
	for (i = 0; i < reg_max; i++)printf("%2d ", odg[i]); printf("\n");
}
void town::get_copy(town* a, int nseq, int reg_max)
{
	int i;
	a->size = size;
	a->fit = fit;
	a->mah = mah;
	a->fpr = fpr;
	for (i = 0; i < size; i++)tot[i].get_copy(&a->tot[i]);
	for (i = 0; i < 16; i++)a->deg[i] = deg[i];
	i = 0;
	for (i = 0; i < reg_max; i++)
	{
		a->odg[i] = odg[i];
	}
	for (i = 0; i < nseq; i++)
	{
		a->pos[i] = pos[i];
		a->ori[i] = ori[i];
	}
}
void town::swap(int n1, int n2)
{
	uno buf;
	tot[n1].get_copy(&buf);
	tot[n2].get_copy(&tot[n1]);
	buf.get_copy(&tot[n2]);
}
int town::order(int n)
{
	/*printf("\n");
	printf("old %d\t",n);
	tot[n].print_all();printf("\n");
	print_all();
	*/
	int ret = n;
	int i;
	int di = 0;
	if (n != 0 && n != size - 1)
	{
		if (tot[n].num < tot[n - 1].num)di = -1;
		else if (tot[n].num > tot[n + 1].num)di = 1;
		else
		{
			if (tot[n].num == tot[n - 1].num && tot[n].sta < tot[n - 1].sta)di = -1;
			if (tot[n].num == tot[n + 1].num && tot[n].sta > tot[n + 1].sta)di = 1;
		}
	}
	{
		if (n == 0)
		{
			if (tot[n].num > tot[n + 1].num)di = 1;
			else if (tot[n].num == tot[n + 1].num && tot[n].sta > tot[n + 1].sta)di = 1;
		}
		if (n == size - 1)
		{
			if (tot[n].num < tot[n - 1].num)di = -1;
			else if (tot[n].num == tot[n - 1].num && tot[n].sta < tot[n - 1].sta)di = -1;
		}
	}
	//printf("DIRECT %d\t",di);
	//if(di==0)return -1;
	if (di == -1)
	{
		for (i = 0; i <= n - 1; i++)
		{
			if (tot[n - i - 1].num > tot[n - i].num)
			{
				//			tot[n-i-1].print_all();
				swap(n - i - 1, n - i);
				ret--;
			}
			else if (tot[n - i - 1].num == tot[n - i].num && tot[n - i - 1].sta > tot[n - i].sta)
			{
				//			tot[n-i-1].print_all();
				swap(n - i - 1, n - i);
				ret--;
			}
			else break;
		}
	}
	if (di == 1)
	{
		for (i = 0; i <= size - n - 2; i++)
		{
			if (tot[n + i + 1].num < tot[n + i].num)
			{
				//			tot[n+i+1].print_all();
				swap(n + i + 1, n + i);
				ret++;
			}
			else if (tot[n + i + 1].num == tot[n + i].num && tot[n + i + 1].sta < tot[n + i].sta)
			{
				//			tot[n+i+1].print_all();
				swap(n + i + 1, n + i);
				ret++;
			}
			else break;
		}
	}
	/*printf("\n");
	printf("new %d\t",ret);
	tot[ret].print_all();printf("\n");
	print_all();
	printf("\n");*/
	return ret;
}
void town::init_rand_hoxa(int nseq, int oln, int rsize, int reg_max, int* len_octa, int* len, int** octa_prowb, int** octa_prows)
{
	int i, j;

	fit = 0;
	fpr = 1;
	size = rsize;
	int oln1 = oln - 1;
	for (j = 0; j < reg_max; j++)odg[j] = 0;
	odg[reg_max] = -1;
	for (i = 0; i < nseq; i++)
	{
		ori[i] = rand() % 2;
		int inx = octa_prows[i][len_octa[i] - 1];
		int r2 = rand() % inx;
		for (j = 0; j < len_octa[i]; j++)
		{
			if (octa_prows[i][j] > r2)
			{
				pos[i] = octa_prowb[i][j];
				break;
			}
		}
	}
	for (i = 0; i < 16; i++)deg[i] = 0;
	i = 0;
	do
	{
		int r = rand() % 16;
		if (deg[r] == oln1)continue;
		deg[r]++;
		i++;
	} while (i < size);
	int t = 0;
	for (i = 0; i < 16; i++)
	{
		if (deg[i] > 0)
		{
			int take_pos[POPSIZE];
			for (j = 0; j < deg[i]; j++)take_pos[j] = 1;
			for (j = deg[i]; j < oln1; j++)take_pos[j] = 0;
			//		printf("DO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			BigMixI(take_pos, oln1);
			//	printf("PO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			for (j = 0; j < oln1; j++)
			{
				if (take_pos[j] == 1)
				{
					tot[t].num = i;
					tot[t].sta = tot[t].end = j;
					t++;
				}
			}
			t -= deg[i];
			int jend = deg[i] - 1;
			for (j = 0; j < deg[i]; j++)
			{
				int rr = rand() % 2;
				int rlen = reg_max - tot[t].end + tot[t].sta;
				if (rr == 0)
				{
					int prev_pos;
					if (j == 0)prev_pos = -1;
					else prev_pos = tot[t - 1].end;
					int spac = tot[t].sta - prev_pos;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].sta -= sh;
					}
				}
				else
				{
					int next_pos;
					if (j != jend)next_pos = tot[t + 1].sta;
					else next_pos = oln1;
					int spac = next_pos - tot[t].end;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].end += sh;
					}
				}
				odg[tot[t].end - tot[t].sta]++;
				t++;
			}
		}
	}
}
/*void town::init_rand(int nseq, int *len, int oln, int rsize, int reg_max)
{
	int i, j;

	fit = 0;
	fpr = 1;
	size = rsize;
	int oln1 = oln - 1;
	for (j = 0; j < reg_max; j++)odg[j] = 0;
	odg[reg_max] = -1;
	//for (i = 0; i < nseq; i++)printf("%d ",len[i]);
	for (i = 0; i < nseq; i++)
	{
		int lenp = len[i] - oln1;
		pos[i] = rand() % lenp;
		ori[i] = rand() % 2;
	}
	for (i = 0; i < 16; i++)deg[i] = 0;
	i = 0;
	do
	{
		int r = rand() % 16;
		if (deg[r] == oln1)continue;
		deg[r]++;
		i++;
	} while (i < size);
	int t = 0;
	for (i = 0; i < 16; i++)
	{
		if (deg[i] > 0)
		{
			int take_pos[POPSIZE];
			for (j = 0; j < deg[i]; j++)take_pos[j] = 1;
			for (j = deg[i]; j < oln1; j++)take_pos[j] = 0;
			//		printf("DO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			BigMixI(take_pos, oln1);
			//	printf("PO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			for (j = 0; j < oln1; j++)
			{
				if (take_pos[j] == 1)
				{
					tot[t].num = i;
					tot[t].sta = tot[t].end = j;
					t++;
				}
			}
			t -= deg[i];
			int jend = deg[i] - 1;
			for (j = 0; j < deg[i]; j++)
			{
				int rr = rand() % 2;
				int rlen = reg_max - tot[t].end + tot[t].sta;
				if (rr == 0)
				{
					int prev_pos;
					if (j == 0)prev_pos = -1;
					else prev_pos = tot[t - 1].end;
					int spac = tot[t].sta - prev_pos;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].sta -= sh;
					}
				}
				else
				{
					int next_pos;
					if (j != jend)next_pos = tot[t + 1].sta;
					else next_pos = oln1;
					int spac = next_pos - tot[t].end;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].end += sh;
					}
				}
				odg[tot[t].end - tot[t].sta]++;
				t++;
			}
		}
	}
}*/
void town::init_rand_part_hoxa(int nseq, int* xporti, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows)
{
	int i, j, k;

	fit = 0;
	fpr = 1;
	//int oln1 = olen - 1;
	for (j = 0; j < nseq; j++)
	{
		i = xporti[j];
		int r = rand() % nind;
		if (r == 0)
		{
			ori[i] = rand() % 2;
			int inx = octa_prows[i][len_octa[i] - 1];
			int r2 = rand() % inx;
			for (k = 0; k < len_octa[i]; k++)
			{
				if (octa_prows[i][k] > r2)
				{
					pos[i] = octa_prowb[i][k];
					break;
				}
			}
		}
	}
}
/*void town::init_rand_part(int nseq, int *len, int oln, int *xporti, int nind)
{
	int i, j;

	fit = 0;
	int oln1 = oln - 1;
	for (j = 0; j < nseq; j++)
	{
		i = xporti[j];
		int r = rand() % nind;
		if (r == 0)
		{
			int lenp = len[i] - oln1;
			pos[i] = rand() % lenp;
			ori[i] = rand() % 2;
		}
	}
}*/
int town::mem_in(int nseq)
{
	pos = new int[nseq];
	if (pos == NULL) return -1;
	//	posc = new int[nseq];
		//if (posc == NULL) return -1;
	ori = new int[nseq];
	if (ori == NULL) return -1;
	return 1;
}
void town::mem_out(void)
{
	delete[] pos;
	//	delete[] posc;
	delete[] ori;
}
int town::sum(int j)
{
	int i;
	int ret = 0;
	if (j == 0) return 0;
	for (i = 0; i < j; i++)ret += deg[i];
	return ret;
}
void town::fprint_all(char* file, char* add)
{
	int i;
	FILE* out;
	char file_out[40];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "\t//%s\n", file);
	fprintf(out, "%d\t%.12f\t", size, fit);
	fprintf(out, "\n");
	for (i = 0; i < size; i++)
	{
		if (i == 0 || (i != 0 && (tot[i].num != tot[i - 1].num)))
		{
			fprintf(out, "%d\t", deg[tot[i].num]);
		}
		fprintf(out, "[%d;%d]%s\t", tot[i].sta, tot[i].end, s[tot[i].num].oli);
		if (tot[i].num != tot[i + 1].num)fprintf(out, "\n");
	}
	fclose(out);
}
void town::fprint_allfi(char* file, char* add, int len, double sd, double c0, double* buf)
{
	int i;
	FILE* out;
	char file_out[40];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "if(strcmp(site,\"%s\")==0)\n", file);// 0 - one window system
	fprintf(out, "{\n\tcity a1 = {");
	fprintf(out, "\n\t//site,size,len,c0,sd\n\t\"%s\",%d,%d,%f,%f,\n\t", file, size, len, c0, sd);
	fprintf(out, "{//buf,sta,end,num\t");
	for (i = 0; odg[i] != -1; i++)fprintf(out, "%d ", odg[i]); fprintf(out, "\n\t");
	for (i = 0; i < size; i++)
	{
		fprintf(out, "{%9f,%d,%d,%d}", buf[i], tot[i].sta, tot[i].end, tot[i].num);
		if (i == size - 1)
		{
			fprintf(out, "}};\t//%s\n", s[tot[i].num].oli);
			break;
		}
		else fprintf(out, ",");
		if (tot[i].num != tot[i + 1].num)fprintf(out, "\t//%s\n\t", s[tot[i].num].oli);
	}
	fprintf(out, "\ta1.get_copy(&sta[win]);\n\twin++;\n}\n");
	fclose(out);
}
void GetWords(int word, int size0, int size, char* w0)
{
	int i, j, zna, che, otn;
	for (i = size0; i < size0 + size; i++)
	{
		che = i - size0;
		zna = size / 4;
		memset(s[i].oli, 0, sizeof(s[i].oli));
		for (j = 0; j < word; j++)
		{
			otn = che / zna;
			s[i].oli[j] = w0[otn];
			che -= otn * zna;
			zna /= 4;
		}
		s[i].oli[word] = '\0';
		//printf("%d\t%s\n",i,s[i].oli);
	}
}
/*
void GetSost2(char *d, stru5 a)
{
int l, i, j, k;
int t[5]={10,11,20,21,22};
int nt=2;
for(i=0;i<a.size;i++)sost[i]=0;
int len=strlen(d), word;
double x, dx;
int left0, right, left;
*/
void MixLong(char* d, double mo[4], int len)
{	// v mo bylo len kratno len zdes
	int i, j;
	int ww[4];
	for (i = 0; i < 4; i++)ww[i] = (int)(mo[i] * len);
	char word[] = "acgt";
	for (i = 0; i < len; i++)
	{
		int r2 = rand() % (len - i);
		int wcur = 0;
		int r;
		for (j = 0; j < 4; j++)
		{
			if (ww[j] == 0)continue;
			wcur += ww[j];
			if (wcur >= r2)
			{
				r = j;
				ww[j]--;
				break;
			}
		}
		d[i] = word[r];
	}
	d[len] = '\0';
}
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
char* TransStr(char* d)
{
	int i, c, lens;
	lens = (int)strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
int CheckStr(char* file, char* d, int n, int print)
{
	int i, len, ret;
	len = (int)strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		int di = (int)d[i];
		if (strchr("atgcATGC", di) != NULL)continue;
		if (strchr("nN", di) != NULL)
		{
			ret = 0; continue;
		}
		if (print == 1)printf("File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
int IdeLet(char c)
{
	int ret;
	switch (c) {
	case 'a': ret = 0; break;
	case 'c': ret = 1; break;
	case 'g': ret = 2; break;
	case 't': ret = 3; break;
	case 'n': ret = -1; break;
	default: ret = -2;
	}
	return(ret);
}
void MinusStr1(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j0 + 1; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr2(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j - 1; k >= j0; k--)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	for (k = 0; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
int BackMat(int size)
{
	int i, j;
	double buf, b[POPSIZE][POPSIZE];
	if (size == 1)
	{
		uw[0][0] = 1 / uw[0][0];
		return 1;
	}
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < size; i++)
		{
			b[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)b[i][i] = 1;
	for (j = 0; j < size - 1; j++)
	{
		for (i = j + 1; i < size; i++)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr1(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				printf("\n1back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//printf("%g(%d)\n",uw[j][j],j);
				return -1;
				//
			}
			uw[i][j] = 0;
		}
	}
	for (j = size - 1; j > 0; j--)
	{
		for (i = j - 1; i >= 0; i--)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr2(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				printf("\n2back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//printf("%g(%d)\n",uw[j][j],j);
				//
				return -1;
			}
			uw[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)
	{
		buf = uw[i][i];
		if (fabs(buf) < 1e-030)
		{
			printf("\nback3 %g %d\t", uw[i][i], i);
			return -1;
		}
		for (j = 0; j < size; j++)
		{
			uw[i][j] /= buf;
			b[i][j] /= buf;
		}
	}
	/*
	for(i=0;i<size;i++)
	{

	for(j=0;j<size;j++)
	{
	printf("%.2e ",uw[i][j]);
	}
	printf("\n");
	}
	printf("\n");
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	printf("%.1e ",b[i][j]);
	}
	printf("\n");
	}
	printf("\n");
	int n;
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	double x=0;
	for(n=0;n<size;n++)x+=uw0[j][n]*b[n][i];
	printf("%.1e ",x);
	}
	printf("\n");
	}
	*/
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			uw[i][j] = b[i][j];
		}
	}
	return 1;
}
/*
void Test(char **peak_real, int *len, int num, int tryx)
{
char d[SEQLEN];
int k;
printf("Test %d\tlen=%d\n",tryx,len[num]);
for(k=0;k<len[num];k++)
{
d[k]=peak_real[num][k];
printf("%c",d[k]);
}
printf("\n");
}
*/
int ComplStr(char* d)
{
	int i, len;
	len = (int)strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a': {d[i] = 't'; break; }
		case 't': {d[i] = 'a'; break; }
		case 'c': {d[i] = 'g'; break; }
		case 'g': {d[i] = 'c'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
int EvalMahControl(town* a, town_ext* best_sel_ext, int nseq, int nseqb, int n_train, int n_cntrl, int* xporti, int* xportj, double* fp_rate, int& n_cntrl_tot, int& n_any_tot, int*** seq, int*** seq_back, int olen, int* len, int* lenb, double** dav, double** dcv, double* qp, qbs* prc)
{
	int k, n, m, o, b, u, z = n_any_tot;
	double av[POPSIZE], buf[POPSIZE];
	double df[POPSIZE];

	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)df[k] = 0;
	for (b = 0; b < n_train; b++)
	{
		m = xporti[b];
		//if (len[m] < olen)continue;		
		int ori = a->ori[m];
		int pos = a->pos[m];
		double fs[POPSIZE];
		for (k = 0; k < a->size; k++)
		{
			int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
			fs[k] = 0;
			for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
			{
				if (a->tot[k].num == seq[ori][m][n + pos])fs[k]++;
			}
			fs[k] /= rlenk;
		}
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			df[k] += fs[k];
		}
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= n_train;
	for (k = 0; k < a->size; k++)df[k] /= n_train;
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)
		{
			uw[k][n] -= df[k] * df[n];
		}
	}
	for (k = 0; k < a->size; k++)
	{
		uw[k][k] += dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= 2;
	for (k = 0; k < a->size; k++)
	{
		av[k] = (df[k] + dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num]) / 2;
		df[k] = df[k] - dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	if (BackMat(a->size) == -1)
	{
		a->fit = 0;
		for (k = 0; k < n_cntrl; k++)fp_rate[n_cntrl_tot + k] = 0.5;
		return 0;
	}
	a->mah = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] = 0;
		for (n = 0; n < a->size; n++)buf[k] += uw[k][n] * df[n];
		a->mah += buf[k] * df[k];
	}
	double c0 = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] /= a->mah;
		c0 -= av[k] * buf[k];
	}
	best_sel_ext->get_copy(c0, buf, a->size);
	for (b = 0; b < n_cntrl; b++)fp_rate[n_cntrl_tot + b] = 0;
	int n_cnt = 0;
	for (b = 0; b < n_cntrl; b++)//control
	{
		u = xportj[b];
		//if (len[u] < olen)continue;		
		int lenp = len[u] - olen + 1;
		double sco_pos = -1000;
		for (m = 0; m < lenp; m++)
		{			
			for (o = 0; o < 2; o++)
			{
				int gom = 1;
				double sco = c0;
				for (k = 0; k < a->size; k++)
				{
					double fs = 0;
					int rlenk = (a->tot[k].end - a->tot[k].sta + 1);					
					for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
					{
						int sym = seq[o][u][n + m];
						if (sym == -1)
						{
							gom = -1;
							break;
						}
						if (a->tot[k].num == sym)fs++;
					}
					if (gom == -1)break;
					if (fs > 0)
					{
						fs /= rlenk;
						sco += buf[k] * fs;
					}
				}
				if (gom == -1)sco = -1;
				else sco = 1 - fabs(1 - sco);
				if (m == 0 && o == 0)sco_pos = sco;
				else
				{
					if (sco > sco_pos)
					{
						sco_pos = sco;
					}
				}
			}
		}
		qp[n_cnt] = prc[z].q = sco_pos;
		prc[z].n = 1;
		n_cnt++;
		z++;
	}
	qsort(qp, n_cntrl, sizeof(double), compare_qq);
	int nseqn = 0;
	for (b = 0; b < nseqb; b++)
	{
		int lenp = lenb[b] - olen + 1;
		nseqn += 2 * lenp;
		double sco_max = -1000;
		for (o = 0; o < 2; o++)
		{
			for (m = 0; m < lenp; m++)
			{
				double sco = c0;
				int gom = 1;
				for (k = 0; k < a->size; k++)
				{
					double fs = 0;
					int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
					for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
					{
						int sym = seq_back[o][b][n + m];
						if (sym == -1)
						{
							gom = -1;
							break;
						}
						if (a->tot[k].num == sym)fs++;
					}
					if (gom == -1)break;
					if (fs > 0)
					{
						fs /= rlenk;
						sco += buf[k] * fs;
					}
				}
				if (gom == -1)sco = -1;
				else sco = 1 - fabs(1 - sco);
				if (sco_max < sco)sco_max = sco;
				if (sco >= qp[0])
				{
					for (k = 0; k < n_cntrl; k++)
					{
						if (sco >= qp[k])
						{
							fp_rate[n_cntrl_tot + k]++;
						}
						else break;
					}
				}
			}
		}
		prc[z].q = sco_max;
		prc[z].n = 0;
		z++;
	}
	for (k = 0; k < n_cntrl; k++)
	{
		int kk = n_cntrl_tot + k;
		//fp_rate[kk] = 0.0001;
		if (fp_rate[kk] > 0)fp_rate[kk] /= nseqn;
		else fp_rate[kk] = 0.5 / (double)nseqn;
	}
	n_cntrl_tot += n_cntrl;
	return 1;
}
double EvalMahFIT(town* a, int n_train, int octa, int* xporti, int*** seq, int olen, double** dav, double** dcv, double** octa_prow, int* len)//double *hoxa_wei, qbs *qps, double **frp, double *octa_rat,
{
	int k, n, m, b;
	double f1[POPSIZE], df[POPSIZE], f2[POPSIZE];
	//	int olen1 = olen - 1;

	//printf("Enter Fit1 %d\n",n_train);
	//exit(1);
	for (k = 0; k < a->size; k++)
	{
		f2[k] = dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)f1[k] = 0;
	//printf("Enter Fit2 %d\n", n_train);
	for (b = 0; b < n_train; b++)
	{
		m = xporti[b];
		//if (len[m] < olen)continue;
		int ori = a->ori[m];
		int pos = a->pos[m];
		double fs[POPSIZE];
		//printf("Seq %d\t",b);
		for (k = 0; k < a->size; k++)
		{
			int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
			fs[k] = 0;
			for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
			{
				int let = seq[ori][m][n + pos];
				if (let == -1)
				{
					a->fpr = 1;
					a->mah = a->fit = 0;
					return 0;
				}
				if (a->tot[k].num == let)fs[k]++;
			}
			fs[k] /= rlenk;
			//frp[b][k] = fs[k];
			//printf("%d %g\t", k, frp[b][k]);
		}
		//printf("\n");		
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			f1[k] += fs[k];
		}
	}
	//exit(1);
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= n_train;
	for (k = 0; k < a->size; k++)f1[k] /= n_train;
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)
		{
			uw[k][n] -= f1[k] * f1[n];
		}
	}
	for (k = 0; k < a->size; k++)
	{
		uw[k][k] += dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= 2;
	//printf("Av1,2\n");
	for (k = 0; k < a->size; k++)
	{
		//av[k] = (f1[k] + f2[k]) / 2;
		df[k] = f1[k] - f2[k];
		//printf("%g %g\t", av[k], df[k]);
	}
	if (BackMat(a->size) == -1) { a->fit = 0; return 0; }
	a->mah = 0;
	for (k = 0; k < a->size; k++)
	{
		double buf = 0;
		for (n = 0; n < a->size; n++)buf += uw[k][n] * df[n];
		a->mah += buf * df[k];
	}
	//exit(1);
	//printf("Mah %f\n", a->mah);
	/*double c0 = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] /= a->mah;
		c0 -= av[k] * buf[k];
		//printf("%d %g\t", k, buf[k]);
	}*/
	//printf("\nc0 %f\n", c0);
	//printf("Scores\n");//frp = NULL;
	//exit(1);
	/*for (b = 0; b < n_train; b++)//train
	{
		//printf("Seq %d\t", b);
		double sco = c0;
		for (k = 0; k < a->size; k++)
		{
			sco += buf[k] * frp[b][k];
		}
		qps[b].q = 1 - fabs(1 - sco);
		qps[b].n = xporti[b];
		//printf("%d %g\n", qps[b].n, qps[b].q);
	}*/
	//printf("\nScoresSort\n");
	//exit(1);
	//qsort(qps, n_train, sizeof(qps[0]), compare_qbs);
	/*for (b = 0; b < n_train; b++)
	{
		m = qps[b].n;
		printf("sco %f peak %d cep %d pos %d\n", qps[b].q,qps[b].n, a->ori[m], a->pos[m]);
	}*/
	double oct = 0;
	int sum = 0;
	//printf("ScoresWei\n");
	//exit(1);	
	//int olen1 = olen - 1;
	for (b = 0; b < n_train; b++)
	{
		//m = qps[b].n;
		m = xporti[b];
		//	oct += octa_prow[a->ori[m]][m][a->pos[m]] * hoxa_wei[b];
		int posdi;
		if (a->ori[m] == 0)posdi = a->pos[m];
		else posdi = len[m] - olen - a->pos[m];
		double wei = octa_prow[m][posdi];
		if ((posdi < 0 || posdi > len[m] - olen) || (wei < -100 || wei > 100))// || posdi > len[m] - olen1)
		{
			printf("ComplHoxaError peak %d ori %d len %d pos %d posdi %d hoxa %g\n", m, a->ori[m], a->pos[m], len[m], posdi, wei);
			exit(1);
		}
		oct += wei;
		sum++;
		//printf("oct %f sco %f peak %d cep %d pos %d pro %f wei %g\n", oct, qps[b].q, m, a->ori[m], a->pos[m], octa_prow[a->ori[m]][m][a->pos[m]], hoxa_wei[b]);
	}
	oct /= sum;
	//printf("Pow\n");
	oct = pow(10, oct);
	a->fpr = oct;
	a->fit = a->mah * a->fpr;
	//exit(1);
	return a->fit;
}
int MutOlig0(town* a)
{
	int cy = 0, r1, i, gom;
	int num1;
	do
	{
		gom = 0;
		r1 = rand() % a->size;
		num1 = rand() % 15;
		if (num1 >= a->tot[r1].num)num1++;
		int i0 = a->sum(num1);
		for (i = i0; i < i0 + a->deg[num1]; i++)
		{
			//if(a->tot[i].num!=num1)continue;
			//if(i==r1)continue;
			if (a->tot[r1].sta >= a->tot[i].sta && a->tot[r1].sta <= a->tot[i].end) { gom = 1; break; }
			if (a->tot[r1].end >= a->tot[i].sta && a->tot[r1].end <= a->tot[i].end) { gom = 1; break; }
			if (a->tot[i].sta >= a->tot[r1].sta && a->tot[i].sta <= a->tot[r1].end) { gom = 1; break; }
			if (a->tot[i].end >= a->tot[r1].sta && a->tot[i].end <= a->tot[r1].end) { gom = 1; break; }
		}
		cy++;
		if (cy > 10)return -1;
	} while (gom == 1);
	if (gom == 1)return -1;
	a->deg[a->tot[r1].num]--;
	a->deg[num1]++;
	a->tot[r1].num = num1;
	return r1;
}
int CTown(town* a, int* n, int len, int olen)
{
	int j, k;
	n[0] = n[1] = -1;
	int i, gom, ns = 0;
	int odgs[] = { 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0 };
	len -= 1;
	for (j = 0; j < a->size; j++)
	{
		if (a->tot[j].sta > a->tot[j].end)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
		if (a->tot[j].end - a->tot[j].sta + 1 > olen)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
	}
	for (j = 0; j < a->size; j++)
	{

		if (a->tot[j].sta >= len || a->tot[j].sta < 0) { n[0] = j; return -1; }
		if (a->tot[j].end >= len || a->tot[j].end < 0) { n[0] = j; return -1; }
		odgs[a->tot[j].end - a->tot[j].sta]++;
	}
	for (i = 0; i < olen; i++)
	{
		if (a->odg[i] != odgs[i])
		{
			printf("Oli %d Announced %d Really %d", i + 1, a->odg[i], odgs[i]);
			return -1;
		}
	}
	for (j = 0; j < 16; j++)
	{
		if (a->deg[j] == 0)continue;
		//if(a->deg[j]==1){ns++;continue;}
		gom = 0;
		for (i = ns; i < ns + a->deg[j]; i++)
		{
			for (k = ns; k < ns + a->deg[j]; k++)
			{
				if (i == k)continue;
				if (a->tot[i].sta >= a->tot[k].sta && a->tot[i].sta <= a->tot[k].end)
				{
					gom = 1; n[0] = i; n[1] = k; break;
				}
				if (a->tot[i].end >= a->tot[k].sta && a->tot[i].end <= a->tot[k].end)
				{
					gom = 1; n[0] = i; n[1] = k; break;
				}
			}
			if (gom == 1)break;
		}
		if (gom == 1)
		{
			printf("Overlap %d\t%d\n", n[0], n[1]);
			a->tot[n[0]].print_all();
			a->tot[n[1]].print_all();
			return -1;
		}
		else ns += a->deg[j];
	}
	int sum = 0;
	for (i = 0; i < olen; i++)
	{
		if (a->odg[i] < 0)
		{
			printf("Oli %d\todg = %d\n", i, a->odg[i]);
			return -1;
		}
		sum += a->odg[i];
	}
	if (sum != a->size)
	{
		printf("Sum by odg %d\tsize = %d\n", sum, a->size);
		return -1;
	}
	return 1;
}
int MutCry0(town* a, int oln, int reg_max)
{
	int n1;
	int a1 = 0, b1 = oln - 1, oln2 = oln - 2;
	int fl = 1, cy = 0;
	do
	{
		n1 = rand() % a->size;
		int n1m1 = n1 - 1;
		if (n1 != 0 && a->tot[n1m1].num == a->tot[n1].num)
		{
			if (a->tot[n1m1].end < oln2)a1 = a->tot[n1m1].end + 1;
			else
			{
				fl = -1;
				cy++;
				continue;
			}
		}
		else a1 = 0;
		int n1p1 = n1 + 1;
		if (n1 != a->size - 1 && a->tot[n1p1].num == a->tot[n1].num)
		{
			if (a->tot[n1p1].sta >= 1)b1 = a->tot[n1p1].sta - 1;
			else
			{
				fl = -1;
				cy++;
				continue;
			}
		}
		else
		{
			b1 = oln2;
		}
		if (cy > 10)return -1;
	} while (fl == -1);
	if (a1 >= b1) return -1;
	int a2, b2, ab = Min(reg_max, b1 - a1 + 1);//max available len oliga 
	int w = rand() % ab;
	int loc = b1 - a1 - w + 1;
	if (loc <= 0)return -1;
	int da2 = rand() % loc;
	a2 = a1 + da2;
	b2 = a2 + w;
	if (a->tot[n1].sta == a2 && a->tot[n1].end == b2)return -1;
	int w0 = a->tot[n1].end - a->tot[n1].sta;
	//	printf("%d%d ",a->odg[w0],a->odg[w]);
	a->odg[w0]--;
	a->odg[w]++;
	//	printf("%d%d \n",a->odg[w0],a->odg[w]);
	a->tot[n1].sta = a2;
	a->tot[n1].end = b2;
	return n1;
}
int MutRegShift(town* a, int n_train, int* len, int* xporti, int olen, int& npeak, int& nori, int& npos)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, cy = 0;
	r1 = rand() % n_train;
	r1 = xporti[r1];
	npeak = r1;
	int lenr = len[r1] - olen + 1;
	r2 = rand() % lenr;
	npos = r2;
	if (r2 == a->pos[r1])
	{
		nori = a->ori[r1] = 1 - a->ori[r1];
		//printf("Out1 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
	else
	{
		a->pos[r1] = r2;
		int r3 = rand() % 2;
		if (r3 == 1)nori = a->ori[r1] = 1 - a->ori[r1];
		else nori = a->ori[r1];
		//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
}
/*int MutRegShiftHoxa(town *a, int n_train, int *xporti, int &npeak, int &nori, int &npos, int *len_octa, int **octa_prowb, int *len, int olen)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, oln1 = olen -1;
	r1 = rand() % n_train;
	npeak = xporti[r1];
	int cep = a->ori[npeak], cep1 = 1 - cep;
	int inx = len_octa[npeak] + len_octa[npeak] - 1;
	r2 = rand() % inx;
	if (r2 < len_octa[npeak])
	{
		nori = cep1;
		npos = octa_prowb[npeak][r2];
	}
	else
	{
		nori = cep;
		r2 -= len_octa[npeak];
		int r3 = octa_prowb[npeak][r2];
		if (r3 >= a->pos[npeak])npos = octa_prowb[npeak][r2 + 1];
		else npos = r3;
	}
	a->ori[npeak] = nori;
	a->pos[npeak] = npos;
	//a->posc[npeak] = len[npeak] - olen - npos;
	//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
	return 1;

}*/
int MutRegShiftHoxaW(town* a, int n_train, int* xporti, int& npeak, int& nori, int& npos, int* len_octa, int** octa_prowb, int** octa_prows)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, j;
	r1 = rand() % n_train;
	npeak = xporti[r1];
	int inx = octa_prows[npeak][len_octa[npeak] - 1];
	r2 = rand() % inx;
	for (j = 0;j < len_octa[npeak];j++)
	{
		if (octa_prows[npeak][j] > r2)
		{
			npos = octa_prowb[npeak][j];
			break;
		}
	}
	if (npos == a->pos[npeak])nori = 1 - a->ori[npeak];
	else nori = rand() % 2;
	a->ori[npeak] = nori;
	a->pos[npeak] = npos;
	//a->posc[npeak] = len[npeak] - olen - npos;
	//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
	return 1;
}
int RecFeat(town a1, town a2, int(*cop)[2], int max)
{
	int ret = 0;
	int reg[2] = { a1.size / 16, a2.size / 16 };
	int cy, r1, r2, num1, num2, i, gom1, gom2, r0;
	do
	{
		cy = 0;
		do
		{
			gom1 = 0;
			r1 = rand() % (a1.size - ret);
			for (i = 0; i < ret; i++)
			{
				if (r1 >= cop[i][0])r1++;
			}
			num1 = a1.tot[r1].num;
			if (a2.deg[num1] == 0)break;
			int i0 = a2.sum(num1);
			for (i = i0; i < i0 + a2.deg[num1]; i++)
			{
				if (gom1 > 1)break;
				if (a2.tot[i].sta > a1.tot[r1].end)break;
				if (a2.tot[i].end < a1.tot[r1].sta)continue;
				if (a2.tot[i].sta >= a1.tot[r1].sta && a2.tot[i].sta <= a1.tot[r1].end) { gom1++; r0 = i; continue; }
				if (a2.tot[i].end >= a1.tot[r1].sta && a2.tot[i].end <= a1.tot[r1].end) { gom1++; r0 = i; continue; }
				if (a1.tot[r1].sta >= a2.tot[i].sta && a1.tot[r1].sta <= a2.tot[i].end) { gom1++; r0 = i; continue; }
				if (a1.tot[r1].end >= a2.tot[i].sta && a1.tot[r1].end <= a2.tot[i].end) { gom1++; r0 = i; continue; }
			}
			if (gom1 == 1)
			{
				if (a1.tot[r1].end == a2.tot[r0].end && a1.tot[r1].sta == a2.tot[r0].sta)
				{
					gom1++;
				}
			}
			if (cy > reg[0])return ret;
			cy++;
		} while (gom1 > 1);
		cy = 0;
		do
		{
			gom2 = 0;
			if (gom1 == 0)r2 = rand() % (a2.size - ret);
			else r2 = r0;
			for (i = 0; i < ret; i++)
			{
				if (r2 >= cop[i][1])r2++;
			}
			num2 = a2.tot[r2].num;
			if (a1.deg[num2] == 0)break;
			if (num1 == num2 && (a1.tot[r1].end == a2.tot[r2].end && a1.tot[r1].sta == a2.tot[r2].sta))continue;
			int i0 = a1.sum(num2);
			for (i = i0; i < i0 + a1.deg[num2]; i++)
			{
				if (gom2 > 1)break;
				if (a1.tot[i].sta > a2.tot[r2].end)break;
				if (a1.tot[i].end < a2.tot[r2].sta)continue;
				if (a1.tot[i].sta >= a2.tot[r2].sta && a1.tot[i].sta <= a2.tot[r2].end) { gom2++; continue; }
				if (a1.tot[i].end >= a2.tot[r2].sta && a1.tot[i].end <= a2.tot[r2].end) { gom2++; continue; }
				if (a2.tot[r2].sta >= a1.tot[i].sta && a2.tot[r2].sta <= a1.tot[i].end) { gom2++; continue; }
				if (a2.tot[r2].end >= a1.tot[i].sta && a2.tot[r2].end <= a1.tot[i].end) { gom2++; continue; }
			}
			if (gom1 == 1 && gom2 == 1)break;
			if (cy > reg[1])return ret;
			cy++;
		} while (gom2 != 0);
		cop[ret][0] = r1;
		cop[ret][1] = r2;
		//		a1.tot[r1].print_all();
		//		a2.tot[r2].print_all();
		//		printf("\n");
		ret++;
	} while (ret < max);
	//	a1.print_all();
	//	a2.print_all();
	return ret;
}
//original: simple exchange
int Reco2_Original(town* a1, town* a2, int* nsi, int* num)
{
	int cy = 0, r1, r2, num1, num2, i, gom1, gom2, r0;
	//int size1=Min(a1->size,a2->size);
	do
	{
		gom1 = 0;
		r1 = rand() % a1->size;
		num1 = a1->tot[r1].num;
		if (a2->deg[num1] == 0)break;
		int i0 = a2->sum(num1);
		for (i = i0; i < i0 + a2->deg[num1]; i++)
		{
			if (gom1 > 1)break;
			if (a2->tot[i].sta > a1->tot[r1].end)break;
			if (a2->tot[i].end < a1->tot[r1].sta)continue;
			if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end) { gom1++; r0 = i; continue; }
			if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end) { gom1++; r0 = i; continue; }
			if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end) { gom1++; r0 = i; continue; }
			if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end) { gom1++; r0 = i; continue; }
		}
		if (gom1 == 1)
		{
			if (a1->tot[r1].end == a2->tot[r0].end && a1->tot[r1].sta == a2->tot[r0].sta)
			{
				gom1++;
			}
		}
		if (cy > 10)return -1;
		cy++;
	} while (gom1 > 1);
	cy = 0;
	do
	{
		gom2 = 0;
		if (gom1 == 0)r2 = rand() % a2->size;
		else r2 = r0;
		num2 = a2->tot[r2].num;
		if (a1->deg[num2] == 0)break;
		if (num1 == num2 && (a1->tot[r1].end == a2->tot[r2].end && a1->tot[r1].sta == a2->tot[r2].sta))continue;
		int i0 = a1->sum(num2);
		for (i = i0; i < i0 + a1->deg[num2]; i++)
		{
			if (gom2 > 1)break;
			if (a1->tot[i].sta > a2->tot[r2].end)break;
			if (a1->tot[i].end < a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end) { gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end) { gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end) { gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end) { gom2++; continue; }
		}
		if (gom1 == 1 && gom2 == 1)break;
		if (cy > 10)return -1;
		cy++;
	} while (gom2 != 0);
	a1->deg[num1]--;
	a2->deg[num2]--;
	a2->deg[num1]++;
	a1->deg[num2]++;
	int w1 = a1->tot[r1].end - a1->tot[r1].sta;
	int w2 = a2->tot[r2].end - a2->tot[r2].sta;
	a1->odg[w1]--;
	a1->odg[w2]++;
	a2->odg[w2]--;
	a2->odg[w1]++;
	uno buf;
	a1->tot[r1].get_copy(&buf);
	a2->tot[r2].get_copy(&a1->tot[r1]);
	buf.get_copy(&a2->tot[r2]);
	num[0] = num1;
	num[1] = num2;
	nsi[0] = r1;
	nsi[1] = r2;
	return 1;
}
//economic: two slightly differnt feature exchanged
int Reco2_Economic(town* a1, town* a2, int* nsi, int* num)
{
	int i0, j, r1, r2, num1, num2, i, gom1 = 0, gom2 = 0, ret;
	int size1 = Min(a1->size, a2->size);
	int ord[POPSIZE];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	size1 /= 2;
	for (j = 0; j < size1; j++)
	{
		r1 = ord[j];
		gom1 = 0;
		num1 = a1->tot[r1].num;
		if (a2->deg[num1] == 0)continue;
		i0 = a2->sum(num1);
		if (a2->deg[num1] == 1)r2 = i0;
		else
		{
			for (i = i0; i < i0 + a2->deg[num1]; i++)
			{
				if (a2->tot[i].end < a1->tot[r1].sta)continue;
				if (gom1 > 1)break;
				if (a2->tot[i].sta > a1->tot[r1].end)break;
				if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end) { gom1++; r2 = i; continue; }//1s 2s 1e
				if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end) { gom1++; r2 = i; continue; }//1s 2e 1e
				if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end) { gom1++; r2 = i; continue; }//2s 1s 2e
				if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end) { gom1++; r2 = i; continue; }//2s 1e 2e
			}
		}
		if (gom1 > 1)continue;
		if (gom1 == 0)
		{
			gom1++;
			r2 = i0 + rand() % a2->deg[num1];
		}
		if (a1->tot[r1].num != a2->tot[r2].num)continue;
		if (a1->tot[r1].end == a2->tot[r2].end && a1->tot[r1].sta == a2->tot[r2].sta)continue;
		/*		if(r2!=0 && a2->tot[r2-1].end>=a1->tot[r1].sta)continue;
		if(r1!=0 && a1->tot[r1-1].end>=a2->tot[r2].sta)continue;
		if(r1!=a1->size-1 && a1->tot[r1+1].sta<=a2->tot[r2].end)continue;
		if(r2!=a2->size-1 && a2->tot[r2+1].sta<=a1->tot[r1].end)continue;
		*/
		num2 = a2->tot[r2].num;
		gom2 = 1;
		i0 = a1->sum(num2);
		for (i = i0; i < i0 + a1->deg[num2]; i++)
		{
			if (i == r1)continue;
			if (gom2 > 1)break;
			if (a1->tot[i].sta > a2->tot[r2].end)break;
			if (a1->tot[i].end < a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end) { gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end) { gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end) { gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end) { gom2++; continue; }
		}
		ret = j;
		if (gom1 == 1 && gom2 == 1)break;
	}
	if (gom1 != 1 || gom2 != 1)return -1;
	if ((r1<0 || r1>a1->size - 1) || (r2<0 || r2>a2->size - 1))
	{
		puts("Error");
	}
	a1->deg[num1]--;
	a2->deg[num2]--;
	a2->deg[num1]++;
	a1->deg[num2]++;
	int w1 = a1->tot[r1].end - a1->tot[r1].sta;
	int w2 = a2->tot[r2].end - a2->tot[r2].sta;
	a1->odg[w1]--;
	a1->odg[w2]++;
	a2->odg[w2]--;
	a2->odg[w1]++;
	uno buf;
	a1->tot[r1].get_copy(&buf);
	a2->tot[r2].get_copy(&a1->tot[r1]);
	buf.get_copy(&a2->tot[r2]);
	num[0] = num1;
	num[1] = num2;
	nsi[0] = r1;
	nsi[1] = r2;
	return ret;
}
//total: exchange of subsets of one dinucl.type
int Reco2_One_dinucleotide_full(town* a1, town* a2)
{
	int gom = -1, i, j, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j < size1; j++)
	{
		if (a1->deg[ord[j]] == 0)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i < a1->deg[ord[j]]; i++)
		{
			int k1 = i1 + i;
			int k2 = i2 + i;
			if (a1->tot[k1].sta != a2->tot[k2].sta) { gom = ord[j]; break; }
			if (a1->tot[k1].end != a2->tot[k2].end) { gom = ord[j]; break; }
		}
		if (gom != -1)break;
	}
	if (gom == -1)return -1;
	uno buf;
	for (j = 0; j < a1->deg[gom]; j++)
	{
		int k1 = i1 + j;
		int k2 = i2 + j;
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a1->odg[w2]++;
		a2->odg[w2]--;
		a2->odg[w1]++;
		a1->tot[k1].get_copy(&buf);
		a2->tot[k2].get_copy(&a1->tot[k1]);
		buf.get_copy(&a2->tot[k2]);
	}
	return 1;
}
//canonical: : classic recombination within one dinucl.type
int Reco2_One_dinucleotide_local(town* a1, town* a2, int olen, int reg_max)
{
	int i, j, j0, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	int rpt = 0, rp[SEQLEN];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j < size1; j++)
	{
		if (a1->deg[ord[j]] <= 1)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i < a1->deg[ord[j]]; i++)
		{
			int k1 = i1 + i;
			int k2 = i2 + i;
			if (a1->tot[k1].sta != a2->tot[k2].sta || a1->tot[k1].end != a2->tot[k2].end)
			{
				if ((a1->tot[k1].sta >= a2->tot[k2].sta && a1->tot[k1].sta <= a2->tot[k2].end) || (a1->tot[k1].end >= a2->tot[k2].sta && a1->tot[k1].end <= a2->tot[k2].end))
				{
					int sta1 = Max(a1->tot[k1].sta, a2->tot[k2].sta);
					int end1 = Min(a1->tot[k1].end, a2->tot[k2].end);
					int over = end1 - sta1;
					if (over < 0)continue;
					int win = Max(a2->tot[k2].end - a1->tot[k1].sta, a1->tot[k1].end - a2->tot[k2].sta);
					if (win >= olen)continue;
					j0 = ord[j];
					rp[rpt++] = i;
				}
			}
		}
		if (rpt != 0)break;
	}
	if (rpt == 0)return -1;
	int rpt0 = rand() % rpt;
	{
		int k1 = i1 + rp[rpt0];
		int k2 = i2 + rp[rpt0];
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a2->odg[w2]--;
		w1 = a2->tot[k2].end - a1->tot[k1].sta;
		if (w1 < 0 || w1 >= reg_max)return -1;
		w2 = a1->tot[k1].end - a2->tot[k2].sta;
		if (w2 < 0 || w2 >= reg_max)return -1;
		a1->odg[w1]++;
		a2->odg[w2]++;
		int end1 = a1->tot[k1].end;
		a1->tot[k1].end = a2->tot[k2].end;
		a2->tot[k2].end = end1;
	}
	uno buf;
	int cha = 0;
	for (j = rp[rpt0] + 1; j < a1->deg[j0]; j++)
	{
		int k1 = i1 + j;
		int k2 = i2 + j;
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a1->odg[w2]++;
		a2->odg[w2]--;
		a2->odg[w1]++;
		a1->tot[k1].get_copy(&buf);
		a2->tot[k2].get_copy(&a1->tot[k1]);
		buf.get_copy(&a2->tot[k2]);
	}
	return 1;
}
int Reco2Peak(town* a1, town* a2, int n_train, int* xporti)
{
	int r1, r2 = 0;
	r1 = rand() % n_train;
	do
	{
		r2 = rand() % n_train;
	} while (r1 == r2);
	r1 = xporti[r1];
	r2 = xporti[r2];
	if (a1->pos[r1] == a2->pos[r2] && a1->ori[r1] == a2->ori[r2])return -1;
	MixI(&a1->ori[r1], &a2->ori[r1]);
	MixI(&a1->pos[r1], &a2->pos[r1]);
	//	MixI(&a1->posc[r1], &a2->posc[r1]);
	MixI(&a1->ori[r2], &a2->ori[r2]);
	MixI(&a1->pos[r2], &a2->pos[r2]);
	//	MixI(&a1->posc[r2], &a2->posc[r2]);
	return 1;
}
void MixPop(town* a, town* b)
{
	town c = *a;
	*a = *b;
	*b = c;
}
int GomTown(town a, town b, int nseq, int* xporti, int check_lpd)
{
	int i, j;
	for (i = 0; i < nseq; i++)
	{
		j = xporti[i];
		if (b.pos[j] != a.pos[j])return 0;
		if (b.ori[j] != a.ori[j])return 0;
	}
	if (check_lpd == 0)return -1;
	for (i = 0; i < a.size; i++)
	{
		if (b.tot[i].sta != a.tot[i].sta)return 0;
		if (b.tot[i].end != a.tot[i].end)return 0;
		if (b.tot[i].num != a.tot[i].num)return 0;
	}
	return -1;
}
int GomTown2(town a, town b)
{
	int i;
	for (i = 0; i < a.size; i++)
	{
		if (b.tot[i].sta != a.tot[i].sta)return 0;
		if (b.tot[i].end != a.tot[i].end)return 0;
		if (b.tot[i].num != a.tot[i].num)return 0;
	}
	return -1;
}
void DelChar(char* str, char c)
{
	int i, lens, size;

	size = 0;
	lens = (int)strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int Fun(char* d, town sta, int len0, double* p, double* rec_buf, double rec_c)
{
	int i, j, ret, len;
	int k;
	char d1[200];

	len = (int)strlen(d);
	ret = 1;
	for (k = 0; k <= len - len0; k++)
	{
		for (j = 0; j < len0; j++)d1[j] = d[k + j]; d1[len0] = '\0';
		p[k] = 0;
		double sco = rec_c;
		for (j = 0; j < sta.size; j++)
		{
			int rlenj = (sta.tot[j].end - sta.tot[j].sta + 1);
			double fm = 0;
			for (i = sta.tot[j].sta; i <= sta.tot[j].end; i++)
			{
				int c1 = IdeLet(d1[i]);
				int c2 = IdeLet(d1[i + 1]);
				if (c1 < 0 || c2 < 0)
				{
					p[k] = 0;
					return 0;
				}
				int cod = 4 * c1 + c2;
				if (sta.tot[j].num == cod) { fm++; break; }
			}
			if (fm != 0)
			{
				fm /= rlenj;
				sco += rec_buf[j] * fm;
			}
		}
		p[k] = 1 - fabs(sco - 1);
	}
	return 1;
}
void GetSost(char* d, int word, int size, int* sost)
{
	int i, j, k, i_sost, let;
	char letter[5] = "acgt";
	int ten[8] = { 1, 4, 16, 64, 256, 1024, 4096, 16384 };
	int lens = (int)strlen(d);
	for (i = 0; i < size; i++)sost[i] = 0;
	int word1 = word - 1;
	int befo = 0, test = 0;
	for (i = 0; i < lens - word1; i++)
	{
		i_sost = 0;
		int coden = IdeLet(d[i + word1]);
		if (befo == 1 && coden >= 0)test = 1;
		if (befo == 0)
		{
			test = 1;
			for (j = 0; j < word; j++)
			{
				int codei = IdeLet(d[i + j]);
				if (codei < 0)
				{
					test = 0;
					break;
				}
			}
		}
		if (test == 1)
		{
			let = 0;
			for (j = word1; j >= 0; j--)
			{
				for (k = 0; k < 4; k++)
				{
					if (d[i + j] == letter[k]) { let = k; break; }
				}
				i_sost += ten[word - 1 - j] * let;
			}
			sost[i_sost]++;
		}
	}
}
void GetSostPro(char* d, int word, int size, int* sost, double* sost_pro)
{
	int i, j, k, i_sost;
	char letter[6] = "acgtn";
	int ten[8] = { 1, 4, 16, 64, 256, 1024, 4096, 16384 };
	int lens = (int)strlen(d);
	int word1 = word - 1;
	for (i = 0; i < size; i++)sost[i] = 0;
	for (i = 0; i < lens - word1; i++)
	{
		i_sost = 0;
		int err = 1;
		for (j = word1; j >= 0; j--)
		{
			for (k = 0; k < 4; k++)
			{
				if (d[i + j] == letter[k])
				{
					i_sost += ten[word1 - j] * k;
					err = 0;
					break;
				}
			}
			if (err == 1)break;
		}
		if (err == 0)
		{
			sost[i_sost]++;
			sost_pro[i] = (double)i_sost;
		}
		else sost_pro[i] = -1;
	}
}
void EvalSeq(char* file, int& nseq, int olen, int len_peak_max)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	nseq = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, n, 1);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)nseq++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
			n++;
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
void EvalLen(char* file, int* len, int olen, int len_peak_max)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalLen! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int nn = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, nn, 0);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)len[n++] = lenx;
			nn++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
void ReadSeq(char* file, int nseq, int* len, int*** seq_real, int olen, double* octaf, int* octa1, int octa, int octa_size, double** octa_pro1, int len_peak_max)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, i, j;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int size = 1;
	for (i = 0; i < octa; i++)size *= 4;

	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d[0]);
			int check = CheckStr(file, d[0], nn, 0);
			nn++;
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)
			{
				TransStr(d[0]);
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				GetSostPro(d[0], octa, octa_size, octa1, octa_pro1[n]);
				for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
				GetSost(d[1], octa, octa_size, octa1);
				for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
				for (j = 0; j < 2; j++)
				{
					int cod[2];
					cod[0] = IdeLet(d[j][0]);
					int len1 = len[n] - 1;
					for (i = 0; i < len1; i++)
					{
						cod[1] = IdeLet(d[j][i + 1]);
						if (cod[0] >= 0 && cod[1] >= 0)seq_real[j][n][i] = 4 * cod[0] + cod[1];
						else
						{
							seq_real[j][n][i] = -1;
						//	printf("ReadSeq! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, j, i, d[j]);
							//exit(1);
						}
						cod[0] = cod[1];
					}
				}
				n++;
			}
			else
			{
				if (lenx < olen)
				{
					printf("ReadSeq! Short peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (lenx > len_peak_max)
				{
					printf("ReadSeq! Long peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (check == -1)
				{
					printf("ReadSeq! Unusual symbol, peak %d ignored\n%s\n", nn + 1, d[0]);
				}
			}
			if (fl == -1)
			{
				fclose(in);
				double sumf = 0;
				for (i = 0; i < octa_size; i++)
				{
					if (octaf[i] == 0)octaf[i] = 1;
					sumf += octaf[i];
				}
				for (i = 0; i < octa_size; i++)octaf[i] /= sumf;
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d[0], 0, sizeof(d[0]));
			DelHole(l);
			strcpy(d[0], l);
			fl = 1; continue;
		}
		if (strlen(d[0]) + strlen(l) > sizeof(d[0]))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d[0], strlen(d[0]));
			exit(1);
		}
		DelHole(l);
		strcat(d[0], l);
	}
}
void ReadSeqBack(char* file, int nseq, int* len, int*** seq_back, int olen, int reg_max, double** dav, double** dcv, double* octaf, int* octa1, int octa, int octa_size, int len_peak_max)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int m, i, j, k, v, fl = 0;
	int n_backr[MOTLEN];
	double freq[16];
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	for (j = 0; j < reg_max; j++)n_backr[j] = 0;
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, nn, 0);
			nn++;
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)
			{
				TransStr(d);
				d[len[n]] = '\0';
				for (j = 0; j < reg_max; j++)n_backr[j] += (len[n] - j + 1);
				int cod[2];
				for (m = 0; m < 2; m++)
				{
					if (m == 1)ComplStr(d);
					GetSost(d, octa, octa_size, octa1);
					for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
					cod[0] = IdeLet(d[0]);
					for (k = 0; k < len[n] - 1; k++)
					{
						cod[1] = IdeLet(d[k + 1]);
						int code = -1;
						if (cod[0] >= 0 && cod[1] >= 0)code = 4 * cod[0] + cod[1];
					//	else
					//	{
						//	printf("ReadSeqBack! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, m, k, d);
							//exit(1);
						//}
						seq_back[m][n][k] = code;
						cod[0] = cod[1];
					}					
					for (j = 1; j <= reg_max; j++)
					{
						int j1 = j - 1;
						for (k = 0; k < len[n] - j; k++)
						{
							int gom = 1;
							for (i = 0; i < 16; i++)freq[i] = 0;
							for (v = 0; v < j; v++)
							{
								int sb = seq_back[m][n][k + v];
								if (sb == -1)
								{
									gom = -1;
									break;
								}
								else freq[sb]++;
							}
							if (gom == 1)
							{
								n_backr[j1]++;
								for (i = 0; i < 16; i++)
								{
									if (freq[i] > 0)
									{
										freq[i] /= j;
										dav[j1][i] += freq[i];
										dcv[j1][i] += freq[i] * freq[i];
									}
								}
							}
						}
					}
				}
				n++;
			}
			else
			{
				if (lenx < olen)
				{
					printf("ReadSeqBack! Short peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (lenx > len_peak_max)
				{
					printf("ReadSeqBack! Long peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (check != 1)
				{
					printf("ReadSeqBack!Unusual symbol, peak %d ignored\n%s\n", nn + 1, d);
					//exit(1);
				}
			}
			if (fl == -1)
			{
				for (j = 0; j < reg_max; j++)
				{
					n_backr[j] *= 2;
					for (k = 0; k < 16; k++)
					{
						dav[j][k] /= n_backr[j];
						dcv[j][k] /= n_backr[j];
						dcv[j][k] -= dav[j][k] * dav[j][k];
					}
				}
				fclose(in);
				double sumf = 0;
				for (i = 0; i < octa_size; i++)
				{
					if (octaf[i] == 0)octaf[i] = 1;
					sumf += octaf[i];
				}
				for (i = 0; i < octa_size; i++)octaf[i] /= sumf;
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
int main(int argc, char* argv[])
{
	int* len, nseq, nseqb, * lenb, i, j, k, n, m;
	char file_for[500], file_back[500], path_fasta[500], path_out[500], pfile_for[500], pfile_back[500];
	int*** seq_real, *** seq_back;
	double** dav;//dinucl.content background
	double** dcv;//self covariations for regions LPD
	//	double **frp;//LPD frequencies
	double* qp;//train scores	
	int** octa_prowb, * len_octa, ** octa_prows;// octa position lists, octa position counts, weight sums
	double** octa_pro1, ** octa_prow, * thr_octa;// , *hoxa_wei;	

	if (argc != 12)
	{
		puts("Sintax: 1path_both_fasta 2file_forground 3file_background 4int max_LPD_length 567int motif_min,max,dif 8double ratio_cnt_of_all(0=jk <0=odd) 9int num_iterations 10 int olig_background 11path_out");//  5<pop_size>
		exit(1);
	}
	//	printf("One ");
	time_t tnow = time(NULL);
	printf("%s", ctime(&tnow));
	//printf("Two ");

	strcpy(path_fasta, argv[1]);
	strcpy(file_for, argv[2]);
	strcpy(file_back, argv[3]);
	strcpy(pfile_for, path_fasta);
	strcpy(pfile_back, path_fasta);
	strcat(pfile_for, file_for);
	strcat(pfile_back, file_back);
	int reg_max = atoi(argv[4]);// max LPD length
	int olen_min = atoi(argv[5]);// dlina motiva
	int olen_max = atoi(argv[6]);// dlina motiva
	int olen_dif = atoi(argv[7]);// dlina motiva	
	double ratio_train_to_control = atof(argv[8]);
	int iteration = atoi(argv[9]);//total no. of jack-knife test			
	int octa = atoi(argv[10]);
	strcpy(path_out, argv[11]);
	double fp2 = 0.001;// FPR threshold for pAUC	
	int len_peak_max = 1000;
	int olen_min0 = 8;
	srand((unsigned)time(NULL));
	dcv = new double* [reg_max];
	if (dcv == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dcv[i] = new double[16];
		if (dcv[i] == NULL) return -1;
	}
	dav = new double* [reg_max];
	if (dav == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dav[i] = new double[16];
		if (dav[i] == NULL) return -1;
	}
	//	int olen1 = olen_max - 1;
	nseq = nseqb = 0;
	//foreground
	EvalSeq(pfile_for, nseq, olen_min, len_peak_max);
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_for, len, olen_min, len_peak_max);
	//background
	EvalSeq(pfile_back, nseqb, olen_min, len_peak_max);
	lenb = new int[nseqb];
	if (lenb == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_back, lenb, olen_min, len_peak_max);
	//
	seq_real = new int** [2];
	if (seq_real == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_real[i] = new int* [nseq];
		for (j = 0; j < nseq; j++)
		{
			seq_real[i][j] = new int[len[j] - 1];
			if (seq_real[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	seq_back = new int** [2];
	if (seq_back == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_back[i] = new int* [nseqb];
		for (j = 0; j < nseqb; j++)
		{
			seq_back[i][j] = new int[lenb[j] - 1];
			if (seq_back[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nseq; j++)for (k = 0; k < len[j] - 1; k++)seq_real[i][j][k] = -1;
		for (j = 0; j < nseqb; j++)for (k = 0; k < lenb[j] - 1; k++)seq_back[i][j][k] = -1;
	}
	for (j = 0; j < reg_max; j++)
	{
		for (i = 0; i < 16; i++)dcv[j][i] = dav[j][i] = 0;
	}
	int octa_size = 1;
	for (i = 0; i < octa; i++)octa_size *= 4;
	octa_pro1 = new double* [nseq];
	if (octa_pro1 == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_pro1[i] = new double[len[i] - octa + 1];
		if (octa_pro1[i] == NULL) return -1;
	}
	octa_prow = new double* [nseq];
	if (octa_prow == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prow[i] = new double[len[i] - octa + 1];
		if (octa_prow[i] == NULL) return -1;
	}
	int len_max = 0;
	for (i = 0; i < nseq; i++)if (len[i] > len_max)len_max = len[i];
	double* octa_pro1p;
	octa_pro1p = new double[len_max];
	if (octa_pro1p == NULL) { puts("Out of memory..."); exit(1); }
	printf("Start1\n");
	double* octa_rat;
	octa_rat = new double[octa_size];
	if (octa_rat == NULL) { puts("Out of memory..."); exit(1); }
	double* octa_av;
	octa_av = new double[octa_size];
	if (octa_av == NULL) { puts("Out of memory..."); exit(1); }
	int* octa1;
	octa1 = new int[octa_size];
	if (octa1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < octa_size; i++)octa_av[i] = 0;
	for (i = 0; i < octa_size; i++)octa1[i] = 0;
	ReadSeq(pfile_for, nseq, len, seq_real, olen_min, octa_av, octa1, octa, octa_size, octa_pro1, len_peak_max);
	for (i = 0; i < octa_size; i++)octa_rat[i] = log10(octa_av[i]);
	for (i = 0; i < octa_size; i++)octa_av[i] = 0;
	for (i = 0; i < octa_size; i++)octa1[i] = 0;
	ReadSeqBack(pfile_back, nseqb, lenb, seq_back, olen_min, reg_max, dav, dcv, octa_av, octa1, octa, octa_size, len_peak_max);
	for (i = 0; i < octa_size; i++)octa_rat[i] -= log10(octa_av[i]);
	/*printf("Octa_rat");
	for (i = 0; i < octa_size; i++)
	{
		if (i % 16 == 0)printf("\n%d", i);
		printf("\t%g", octa_rat[i]);
	}
	printf("\n");*/
	for (i = 0; i < nseq; i++)
	{
		int leni = len[i] - octa + 1;
		for (n = 0; n < leni; n++)
		{
			int i_sost = (int)octa_pro1[i][n];
			octa_pro1[i][n] = octa_rat[i_sost];
		}
	}
	//exit(1);
	for (j = 0; j < reg_max; j++)
	{
		printf("AvR%d\t", j + 1);
		for (i = 0; i < 16; i++)printf("%.4f\t", dav[j][i]);
		printf("\n");
	}
	for (j = 0; j < reg_max; j++)
	{
		printf("CovR%d\t", j + 1);
		for (i = 0; i < 16; i++)printf("%.4f\t", dcv[j][i]);
		printf("\n");
	}
	//Test(peak_real[0],len,0,1);
	if (ratio_train_to_control == 0)
	{
		iteration = Min(nseq, iteration);
	}
	//	printf("Start 5 Iter %d Nseq %d\n",iteration, nseq);
	int* n_train;
	n_train = new int[iteration];
	if (n_train == NULL) { puts("Out of memory..."); exit(1); }
	int* n_cntrl;
	n_cntrl = new int[iteration];
	if (n_cntrl == NULL) { puts("Out of memory..."); exit(1); }
	int iter;
	if (ratio_train_to_control == 0) //Leave-one-out cross-validation
	{
		int nseq1 = nseq - 1;
		for (i = 0; i < iteration; i++) { n_train[i] = nseq1; n_cntrl[i] = 1; }
	}
	else
	{
		if (ratio_train_to_control > 0)// random splits cross-validation
		{
			int n_cnt = (int)(nseq / ratio_train_to_control);
			int n_trn = nseq - n_cnt;
			for (iter = 0; iter < iteration; iter++)
			{
				n_train[iter] = n_trn;
				n_cntrl[iter] = n_cnt;
			}
		}
		else //two-fold cross-validation
		{
			int n_less = nseq / 2;
			int n_more = nseq - n_less;
			for (iter = 0; iter < iteration; iter++)
			{
				if (iter % 2 == 0)
				{
					n_train[iter] = n_less;
					n_cntrl[iter] = n_more;
				}
				else
				{
					n_train[iter] = n_more;
					n_cntrl[iter] = n_less;
				}
			}
		}
	}
	int n_cnt_tot = 0;
	for (i = 0; i < iteration; i++)n_cnt_tot += n_cntrl[i];
	int n_trn_tot = 0;
	for (i = 0; i < iteration; i++)n_trn_tot += n_train[i];
	int n_both_sam = 0;
	for (i = 0; i < iteration; i++)
	{
		n_both_sam += n_cntrl[i];
		n_both_sam += nseqb;
	}
	qbs* prc;
	prc = new qbs[n_both_sam];
	if (prc == NULL)
	{
		puts("Out of memory..."); exit(1);
	}
	for (i = 0; i < n_both_sam; i++)
	{
		prc[i].n = 0;
		prc[i].q = 0;
	}
	qbs* prc_best;
	prc_best = new qbs[n_both_sam];
	if (prc == NULL)
	{
		puts("Out of memory..."); exit(1);
	}
	for (i = 0; i < n_both_sam; i++)
	{
		prc_best[i].n = 0;
		prc_best[i].q = 0;
	}
	double* fp_rate;
	fp_rate = new double[n_cnt_tot + 10];
	if (fp_rate == NULL)
	{
		puts("Out of memory..."); exit(1);
	}
	int* xport;
	xport = new int[nseq];
	if (xport == NULL) { puts("Out of memory..."); exit(1); }
	int* xporti;//train
	xporti = new int[nseq];
	if (xporti == NULL) { puts("Out of memory..."); exit(1); }
	int* xportj;//control
	xportj = new int[nseq];
	if (xportj == NULL) { puts("Out of memory..."); exit(1); }
	{
		char word[] = "acgt";
		GetWords(2, 0, 16, word);
	}
	town_ext pop_ext;
	pop = new town * [iteration];
	if (pop == NULL) { puts("Out of memory...!"); exit(1); }
	for (iter = 0; iter < iteration; iter++)
	{
		pop[iter] = new town[MEGE];
		if (pop[iter] == NULL) { puts("Out of memory...!"); exit(1); }
	}
	for (iter = 0; iter < iteration; iter++)
	{
		for (i = 0; i < MEGE; i++)
		{
			pop[iter][i].mem_in(nseq);
		}
	}
	for (i = 0; i < 2; i++)
	{
		det2[i].mem_in(nseq);
	}
	det1.mem_in(nseq);
	int n_decil[CENT];
	{
		n_decil[0] = 1;
		for (n = 1; n < CENT; n++)n_decil[n] = (int)(n_cnt_tot * n / CENT);
	}
	double auc_roc_max = 0, auc_pr_max = 0;
	int size_selected, olen, olen_selected = olen_min;
	int size_selected2, olen_selected2 = olen_min;
	double* fp_rate_best;
	fp_rate_best = new double[n_cnt_tot + 10];
	if (fp_rate_best == NULL) { puts("Out of memory..."); exit(1); }
	/*	int *tp_rate;
	tp_rate = new int[n_cnt_tot + 10];
	if (tp_rate == NULL) { puts("Out of memory..."); exit(1); }
	double* fp_rate_step;
	fp_rate_step = new double[n_cnt_tot + 10];
	if (fp_rate_step == NULL) { puts("Out of memory..."); exit(1); }*/
	//	double dtp = 1 / (double)nseq;
	char add_roc[500], add_auc[500], add_prc[500], add_roc1[500], add_prc1[500];
	strcpy(add_roc, "_roc_bs.txt");
	strcpy(add_auc, "_auc_bs.txt");
	strcpy(add_prc1, "_prc_bs1.txt");
	strcpy(add_roc1, "_roc_bs1.txt");
	strcpy(add_prc, "_prc_bs.txt");
	char file_out_cnt[500];
	int n_train_max = 0;
	int size_start = 40;// (int)(k_size_start*olenf_max);
	int size_end = 100;// (int)(k_size_end*olenf_max);
	int size_dif = 20;// (int)(k_size_dif*olenf_max);
	for (iter = 0; iter < iteration; iter++)
	{
		if (n_train[iter] > n_train_max)n_train_max = n_train[iter];
	}
	qp = new double[nseq];
	if (qp == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < nseq; i++)qp[i] = 0;
	thr_octa = new double[nseq];
	if (thr_octa == NULL) { puts("Out of memory..."); exit(1); }
	len_octa = new int[nseq];
	if (len_octa == NULL) { puts("Out of memory..."); exit(1); }
	octa_prowb = new int* [nseq];
	if (octa_prowb == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prowb[i] = new int[len[i]];
		//octa_prowb[i] = new int[len[i]];
		if (octa_prowb[i] == NULL) return -1;
	}
	for (i = 0; i < nseq; i++)for (n = 0; n < len[i]; n++)octa_prowb[i][n] = -1;
	octa_prows = new int* [nseq];
	if (octa_prows == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prows[i] = new int[len[i]];
		//octa_prows[i] = new int[len[i]];
		if (octa_prows[i] == NULL) return -1;
	}
	//printf("Octa_prowb2\n");	
//	tnow = time(NULL);
//	printf("%s", ctime(&tnow));
//	exit(1);
	for (olen = olen_min; olen <= olen_max; olen += olen_dif)
	{
		//printf("Octa_prow Ini1\n");
		int len_tot = 0, len_wei = 0;
		for (i = 0; i < nseq; i++)
		{
			int leni = len[i] - olen + 1;
			len_tot += leni;
			int odif = olen - octa + 1;
			for (n = 0; n < leni; n++)
			{
				octa_prow[i][n] = 0;
				for (k = 0; k < odif; k++)octa_prow[i][n] += octa_pro1[i][n + k];
				octa_prow[i][n] /= odif;
			}
			int half = leni / 3 - 1;
			{
				for (n = 0; n < leni; n++)octa_pro1p[n] = octa_prow[i][n];
				qsort(octa_pro1p, leni, sizeof(double), compare_qq2);
				//if(octa_pro1p[half] < 0)
				thr_octa[i] = octa_pro1p[half];
				//else thr_octa[i] = 0;
			}
			//	double sumw = 0;
			double maxw = 0;
			k = 0;
			for (n = 0; n < leni; n++)
			{
				double dw = octa_prow[i][n] - thr_octa[i];
				if (dw >= 0)
				{
					//sumw += dw;					
					octa_prowb[i][k] = n;
					k++;
					if (dw > maxw)dw = maxw;
				}
			}
			len_octa[i] = k;
			len_wei += len_octa[i];
			double koef;
			if (maxw > 0) { koef = Max(1, 5 / maxw); }
			else koef = 1;
			octa_prows[i][0] = 1 + (int)(koef * (octa_prow[i][octa_prowb[i][0]] - thr_octa[i]));
			for (n = 1; n < len_octa[i]; n++)
			{
				octa_prows[i][n] = octa_prows[i][n - 1] + 1 + (int)(koef * (octa_prow[i][octa_prowb[i][n]] - thr_octa[i]));
			}
			//printf("%d %d\t", i,len_octa[i]);
			//if (i % 20 == 0)printf("\n");
		}
		printf("FractHoxa %f\n", (double)len_wei / len_tot);
		int size0;
		int size_len = size_start;
		int size_len2 = size_start;
		double auc_roc_len = 0, auc_pr_len = 0;
		for (size0 = size_start; size0 <= size_end; size0 += size_dif)
		{
			for (k = 0; k < n_cnt_tot + 1; k++)fp_rate[k] = 0;
			int big_exit1 = 1;// local exit (separ +-) global exit (separation do not exceeded the previous run)
			double fit_prev, fit_after_mut;
			int cnt_count = 0;
			int n_any_tot = 0;
			//Test(peak_real[0],len,0,2);	
			for (iter = 0; iter < iteration; iter++)
			{
				printf("\n%s vs %s\tWindowlength %d\tBootStrap %d  the last %d\t", file_for, file_back, olen, iter + 1, iteration);
				printf("N LPD %d\tDeg %d\tEli %d\n", size0, MEGE, ELIT);
				if (ratio_train_to_control > 0)// random splits cross-validation
				{
					for (k = 0; k < n_cntrl[iter]; k++)xport[k] = 0;
					for (k = n_cntrl[iter]; k < nseq; k++)xport[k] = 1;
					BigMixI(xport, nseq);
				}
				else
				{
					if (ratio_train_to_control == 0) //Leave-one-out cross-validation
					{
						for (k = 0; k < nseq; k++)xport[k] = 1;//train
						xport[iter] = 0;//control
					}
					else //two-fold cross-validation
					{

						if (iter % 2 == 0)
						{
							xport[0] = 0;//control
						}
						else
						{
							xport[0] = 1;//train
						}
						for (k = 1; k < nseq; k++)
						{
							xport[k] = 1 - xport[k - 1];
						}
					}
				}
				{
					i = 0;//train
					for (k = 0; k < nseq; k++)
					{
						if (xport[k] == 1)xporti[i++] = k;
					}
				}
				//printf("Xport Ini2\n");
				//exit(1);
				//initiation		
				int pair_all;
				int pair_d[10000][2];
				int pair_take[20000];
				int gen = 0;
				int success_r, success_m;
				int success_r1[NREC];
				//PARAMETERS SETTING																				
				int check_lpd = 1;
				//	int pair_all1 = pair_all - 1;
				int stop_oi[MEGE], stop_li[MEGE], stop_pi[MEGE];
				int success_ri[MEGE];
				for (i = 0; i < MEGE; i++)stop_oi[i] = stop_li[i] = stop_pi[i] = 0;
				for (i = 0; i < MEGE; i++)success_ri[i] = 1;
				int mege_h;
				int restart = 0;
				//Test(peak_real[0],len,0,3);
				do
				{
					int success_o, success_l, success_p;
					if (big_exit1 == 1)
					{
						for (i = 0; i < MEGE; i++)
						{
							int err = 1;
							//			printf("Poo Ini1 %d\n", i);							
							do
							{
								err = 1;
								int gom = 0;
								det1.init_rand_hoxa(nseq, olen, size0, reg_max, len_octa, len, octa_prowb, octa_prows);
								if (det1.check(0, reg_max) == -1)
								{
									det1.check(0, reg_max);
									printf("Population error!\n");
									det1.print_all(reg_max, nseq);
									exit(1);
								}
								else
								{
									for (m = 0; m < i; m++)
									{
										gom = GomTown(det1, pop[iter][m], n_train[iter], xporti, 1);
										if (gom == -1)
										{
											break;
										}
									}
								}
								if (gom != -1)
								{
									//							printf("Poo Ini2 %d\n", i);
									EvalMahFIT(&det1, n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len);//hoxa_wei,qps, frp, octa_rat,
									det1.print_sta(reg_max);
									det1.get_copy(&pop[iter][i], nseq, reg_max);
									//						printf("Poo Ini3 %d\n", i);
									err = 0;
								}
							} while (err == 1);
						}
						big_exit1 = 0;
						qsort((void*)(&pop[iter][0]), MEGE, sizeof(pop[iter][0]), compare_pop);
						//pop[iter][0].print_all(reg_max,nseq);
					}
					/*	for(i=0;i<MEGE;i++)
					{
					//printf("After %d\n",i+1);
					//	pop[i].print_all(reg_max,nseq);
					if(pop[iter][i].check(0,reg_max)==-1)
					{
					pop[iter][i].print_all(reg_max,nseq);
					printf("Population error!\n");
					exit(1);
					}
					}*/
					//	printf("Poo Ini2\n");
						//exit(1);
					fit_prev = pop[iter][0].fit;
					success_o = success_l = success_p = success_m = 0;
					double ratio2_gen0 = 0.01, ratio_rec_cycle = 0.005, ratio_mut_cycle = 0.0001;
					int step, step_max, step_max_tot = 0;
					int elit_rec;
					int kn_train_max = Max(n_train_max, 500);
					int sr_step = 200 * kn_train_max;//100000 if n_train_max = 500		;
					int n_rec_cycle_max = 2000 * kn_train_max;//1000000 if n_train_max = 500		
					double jwei;
					if (restart == 0)
					{
						step = 2000;
						step_max = 20 * kn_train_max; //10000 if n_train_max = 500						
						elit_rec = ELIT / 4;
						jwei = 1.1;
						mege_h = MEGE;
					}
					else
					{
						if (restart == 1)
						{
							jwei = 1.2;
							step = 25000;
							step_max = 500 * kn_train_max;//250000 if n_train_max = 500							
							elit_rec = ELIT / 5;
							mege_h = ELIT;
						}
						else
						{
							step = 50000;
							step_max = 1000 * kn_train_max;//500000 if n_train_max = 500							
							elit_rec = 0;
							jwei = 1.3;
						}
					}
					if (step_max < step)step_max = step;
					double ratio_thr = 1 / (double)step;
					double ratio_thr_r[2];
					for (i = 0; i < 2; i++)ratio_thr_r[i] = ratio_thr;
					int step2 = 2 * step;
					int mege_h1 = mege_h - 1;
					//mutations									
					int n_mut_tot = 0, success_m_tot = 0, mdo = 1;
					int asuccess[NMUT], atry[NMUT];
					for (k = 0; k < NMUT; k++)asuccess[k] = atry[k] = 0;
					int success_mi[MEGE], try_mi[MEGE];
					double exp_rec_rate[MEGE];
					for (i = 0; i < mege_h; i++)exp_rec_rate[i] = 0;
					for (i = 0; i < mege_h; i++)success_mi[i] = try_mi[i] = 0;
					double sm0_rate1[2] = { 0,0 };
					int mut_jump = 0;
					do
					{
						step_max_tot += step_max;
						asuccess[2] = atry[2] = 0;
						for (i = 0; i < mege_h; i++)pop[iter][i].sm = pop[iter][i].tm = 0;
						if (gen == 0)printf("Mut cycle %d\n", step_max_tot / step_max);
						//double fit_pre_mut = pop[iter][0].fit;						
						for (i = 0; i < mege_h; i++)
						{
							if (stop_pi[i] == 0 || (stop_li[i] == 0 || stop_oi[i] == 0))
							{
								//	printf("Mut start %d\n",i+1);
								pop[iter][i].get_copy(&det1, nseq, reg_max);
								int success_o_local = 0, success_l_local = 0, success_p_local = 0;
								int success_m_local = 0;
								double ratio[NMUT];
								if (stop_oi[i] == 0)ratio[0] = 1;
								else ratio[0] = 0;
								if (stop_li[i] == 0)ratio[1] = 1;
								else ratio[1] = 0;
								if (stop_pi[i] == 0)ratio[2] = 1;
								else ratio[2] = 0;
								int step_success[NMUT];
								int step_try[NMUT];
								for (k = 0; k < NMUT; k++)step_success[k] = step_try[k] = 0;
								int n_mut_here = 0;
								double fit_mut_prev0 = pop[iter][i].fit;
								int n_mut_min_cyc = step;
								int m_iter = 0;
								while (ratio[2] > ratio_thr || (ratio[0] > ratio_thr || ratio[1] > ratio_thr))
								{
									//mm++;
									//printf("I=%d Mut %d\t%f %f %f\t",i+1,mm,ratio[0],ratio[1],ratio[2]);
									//Test(peak_real[0],len,0,4);
									int sm;
									//if (stop_li[i] == 1 && stop_oi[i] == 1)sm = 2;
									//else
									{
										int rs = step + step_success[2];
										if (ratio[1] > ratio_thr && stop_li[i] == 0)rs += step + step_success[1];
										if (ratio[0] > ratio_thr && stop_oi[i] == 0)rs += step + step_success[0];
										int rs1 = rand() % rs;
										if (rs1 < step + step_success[2])sm = 2;
										else
										{
											if (rs1 < step2 + step_success[2] + step_success[1])sm = 1;
											else sm = 0;
										}
										//sm = rand() % 3;
									}
									//if(sm==2 && stop_p==1)continue;
									//if(sm==1 && stop_l==1)continue;
									//if(sm==0 && stop_o==1)continue;																
									int muto, muto1;
									int npeak = -1, nori = -1, npos = -1;
									//printf("Sm=%d Wei %d\n",sm,wei);
									//Test(peak_real[0],len,0,5);					
									if (sm == 0)muto = MutOlig0(&det1);
									else
									{
										if (sm == 1)muto = MutCry0(&det1, olen, reg_max);
										//else muto = MutRegShift(&det1, n_train[iter], len, xporti, olen, npeak, nori, npos);
										//else muto = MutRegShiftHoxa0(&det1, n_train[iter],len,xporti, olen,npeak, nori, npos,thr_octa,octa_prow);
										//else muto = MutRegShiftHoxa(&det1, n_train[iter], xporti, npeak, nori, npos, len_octa, octa_prowb, len, olen);
										else muto = MutRegShiftHoxaW(&det1, n_train[iter], xporti, npeak, nori, npos, len_octa, octa_prowb, octa_prows);
									}
									int gom = 0;
									//Test(peak_real[0],len,0,6);								
									if (muto != -1)
									{
										muto1 = det1.order(muto);
										if (muto1 == -1)
										{
											puts("MUTO Error");
											gom = -1;
										}
										if (gom == 0)
										{
											for (m = 0; m < mege_h; m++)
											{
												if (m == i)continue;
												gom = GomTown(det1, pop[iter][m], n_train[iter], xporti, check_lpd);
												if (gom == -1) { break; }
											}
										}
									}
									//	Test(peak_real[0],len,0,7);					
									//	if(det1.check(0,reg_max)==-1){printf("Population error!\n");exit(1);}
									//det.print_all();				
									double dd = 0;
									if (gom != -1)
									{
										n_mut_here++;
										EvalMahFIT(&det1, n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len);//hoxa_wei,qps, frp, octa_rat,
										dd = det1.fit / pop[iter][i].fit;
										pop[iter][i].tm++;
										if (dd > 1)
										{
											step_success[sm]++;
											pop[iter][i].sm++;
											if (sm == 0)success_o_local++;
											else
											{
												if (sm == 1)success_l_local++;
												else success_p_local++;
											}
											success_m_local++;
											det1.get_copy(&pop[iter][i], nseq, reg_max);
											success_mi[i]++;
										}
									}
									if (dd <= 1)
									{
										pop[iter][i].get_copy(&det1, nseq, reg_max);
									}
									step_try[sm]++;
									try_mi[i]++;
									if (n_mut_here % step == 0)
									{
										m_iter++;
										for (k = 0; k < 3; k++)
										{
											if (step_success[k] != 0 && step_try[k] != 0)
											{
												ratio[k] = (double)step_success[k] / step_try[k];
											}
											else ratio[k] = 0;
										}
										if (ratio[0] <= ratio_thr)stop_oi[i] = 1;
										if (ratio[1] <= ratio_thr)stop_li[i] = 1;
										if (ratio[2] <= ratio_thr)stop_pi[i] = 1;
										double ratio_per_cycle = pop[iter][i].fit / fit_mut_prev0 - 1;
										for (k = 0; k < 3; k++)atry[k] += step_try[k];
										for (k = 0; k < 3; k++)asuccess[k] += step_success[k];
										int step_success_tot = step_success[0] + step_success[1] + step_success[2];
										if (step_success_tot < n_mut_min_cyc)n_mut_min_cyc = step_success_tot;
										if (try_mi[i] != 0)exp_rec_rate[i] = (double)success_mi[i] / try_mi[i];
										else exp_rec_rate[i] = 0;
										if (gen == 1)printf("Step%d M%d %d,%d,%d Try %d,%d,%d Min %d M %f H %g Fit %f Ratio %f RatioSco %f\n", m_iter, i + 1, step_success[0], step_success[1], step_success[2], step_try[0], step_try[1], step_try[2], n_mut_min_cyc, pop[iter][i].mah, pop[iter][i].fpr, pop[iter][i].fit, ratio[2], ratio_per_cycle);
										if (gen > 1)printf("Step%d M%d %d Try %d Min %d M %f H %g Fit %f RatioM %f RatioSco %f\n", m_iter, i + 1, step_success[2], step_try[2], n_mut_min_cyc, pop[iter][i].mah, pop[iter][i].fpr, pop[iter][i].fit, ratio[2], ratio_per_cycle);
										success_mi[i] = try_mi[i] = 0;
										for (k = 0; k < 3; k++)step_try[k] = step_success[k] = 0;
										if (n_mut_here >= step_max)break;
										if (ratio_per_cycle <= ratio_mut_cycle)
										{
											//printf("Too small score growth %f\n", ratio_per_cycle);
											stop_pi[i] = 1;
											break;
										}
										fit_mut_prev0 = pop[iter][i].fit;
									}
								}
								success_o += success_o_local;
								success_l += success_l_local;
								success_p += success_p_local;
								success_m_tot += success_m_local;
								printf("M %d %d,%d,%d = %d Try %d M %f H %g F %f", i + 1, success_o_local, success_l_local, success_p_local, success_m_local, n_mut_here, pop[iter][i].mah, pop[iter][i].fpr, pop[iter][i].fit);
								if (gen == 0)printf("\tOLP %d%d%d RatioP %f", stop_oi[i], stop_li[i], stop_pi[i], ratio[2]);
								printf("\n");
								n_mut_tot += n_mut_here;
								if (success_m_local != 0)success_m++;
							}
						}
						mut_jump = 0;
						for (i = 1; i < mege_h; i++)
						{
							if (pop[iter][i].fit >= pop[iter][0].fit)mut_jump++;
						}
						qsort((void*)(&pop[iter][0]), mege_h, sizeof(pop[iter][0]), compare_pop);
						sm0_rate1[0] = sm0_rate1[1];
						{
							int is = 0, trym = 0, sucm = 0;
							for (i = 0; i < mege_h; i++)
							{
								if (pop[iter][i].tr > 0) { is = i; break; }
							}
							int ie;
							if (gen == 0) { ie = Min(is + elit_rec, mege_h1); }
							else { ie = Min(is + mut_jump, mege_h1); }
							//for (i = is; i <= ie; i++)printf("%d TM %d SM %d\n", i + 1, pop[iter][i].tm, pop[iter][i].sm);
							for (i = is; i <= ie; i++)
							{
								sucm += pop[iter][i].sm;
								trym += pop[iter][i].tm;
							}
							sm0_rate1[1] = (double)sucm / trym;
						}
						if (gen > 0)
						{
							mdo = 0;
							break;
						}
						else
						{
							mdo = 0;
							if (mut_jump >= ELIT)
							{
								mdo = 1;
								printf("Mut jump %d\tSM rate %g\n", mut_jump, sm0_rate1[1]);
							}
							else
							{
								int m_nostoplo = 0, m_stopp = 0;
								for (i = 0; i < mege_h; i++)
								{
									if ((stop_li[i] == 0 || stop_oi[i] == 0) && stop_pi[i] == 0)
									{
										m_nostoplo++;
									}
									if (stop_pi[i] == 1)
									{
										m_stopp++;
									}
								}
								double rat0 = 0;
								if (atry[2] > 0)
								{
									rat0 = (double)asuccess[2] / atry[2];
								}
								else mdo = 0;
								if (m_nostoplo > ELIT || rat0 > ratio2_gen0)mdo = 1;
								if (m_stopp >= ELIT)mdo = 0;
								printf("Common RatioP %f\tMut jump %d\tyLOnP %d yP %d\tSM rate %g\t", rat0, mut_jump, m_nostoplo, m_stopp, sm0_rate1[1]);
								tnow = time(NULL);
								printf("%s", ctime(&tnow));
							}
						}
					} while (mdo == 1);
					if (gen == 0)success_m /= (step_max_tot / step_max);
					if (gen <= 1)ratio_thr_r[0] = ratio_thr;		//OL								
					else ratio_thr_r[0] = 1;
					//	elit_rec = Max(elit_rec, mut_jump);
					//	elit_rec = Min(elit_rec, mege_h / 2 - 1);
					double sm0_rate;
					{
						if (gen == 0)
						{
							sm0_rate = sm0_rate1[1] * sm0_rate1[1] / sm0_rate1[0];
						}
						else sm0_rate = sm0_rate1[1];
					}
					//sm0_rate /= 2;
					//if (gen == 0)sm0_rate /= 4;
					//else  sm0_rate /= 2;
					for (i = 0; i < mege_h; i++)pop[iter][i].sr = pop[iter][i].tr = 0;
					{//P
						int recount = 0;
						double resum = 0;
						for (i = 0; i < mege_h - 1; i++)
						{
							double e1 = Max(ratio_thr, exp_rec_rate[i]);
							double e2 = Max(ratio_thr, exp_rec_rate[i + 1]);
							resum += success_ri[i] * (e1 + e2) / 2;
							recount += success_ri[i];
						}
						if (recount != 0)ratio_thr_r[1] = resum / recount;
						else ratio_thr_r[1] = ratio_thr;
					}
					//pop[iter][0].print_all(reg_max, nseq);
					//recombinations											
					printf("Total mut %d\tSM rate %g\t Mut jump %d\n", n_mut_tot, sm0_rate, mut_jump);
					int loc_rec, loc_rec_tot = 0;
					fit_after_mut = pop[iter][0].fit;
					for (m = 0; m < NREC; m++)success_r1[m] = 0;
					int rwei[MEGE];
					double jwei1 = jwei;
					{
						int jmax = mege_h / 2;
						rwei[jmax] = 2;
						for (j = jmax - 1; j >= 0; j--)
						{
							rwei[j] = (int)(rwei[jmax] * jwei1);
							jwei1 *= jwei;
						}
						for (j = jmax; j >= 0; j--)printf("%d %d\t", j, rwei[j]);
						printf("\n");
						pair_all = 0;
						{
							for (j = jmax; j >= 0; j--)
							{
								for (k = j + 1; k < mege_h; k++)
								{
									for (m = 0; m < rwei[j]; m++)
									{
										pair_d[pair_all][0] = j;
										pair_d[pair_all][1] = k;
										pair_all++;
										pair_d[pair_all][0] = k;
										pair_d[pair_all][1] = j;
										pair_all++;
									}
								}
							}
						}
					}
					for (k = 0; k < pair_all; k++)pair_take[k] = k;
					int n_rec_cycle, sr_step_cycle = sr_step;
					if (restart <= 1) { n_rec_cycle = Min(n_mut_tot / 4, n_rec_cycle_max); }
					else n_rec_cycle = n_rec_cycle_max;
					n_rec_cycle /= pair_all;
					if (n_rec_cycle < 1)n_rec_cycle = 1;
					sr_step_cycle /= pair_all;
					int rec_jump[MEGE];
					for (k = 0; k < mege_h; k++)rec_jump[k] = 0;
					double fit_pre_rec[MEGE];
					for (k = 0; k < mege_h; k++)fit_pre_rec[k] = pop[iter][k].fit;
					int sr;
					double ratio_r[2] = { 1,1 };
					int step_rtry[NREC], step_rsuccess[NREC];
					for (k = 0; k < NREC; k++)step_rtry[k] = step_rsuccess[k] = 0;
					int success_r_cycle = 0;
					printf("Rec %d cycles of %d tries\tTotal %d\tRatioThrOL %.5f RatioThrP %.5f StepR %d\t", n_rec_cycle, pair_all, n_rec_cycle * pair_all, ratio_thr_r[0], ratio_thr_r[1], sr_step_cycle);
					printf("RE %d", elit_rec);
					printf("\n");
					double fit_rec_prev = pop[iter][elit_rec].fit;
					double fit_rec_prev0 = pop[iter][0].fit;
					double mah_rec_prev = pop[iter][0].mah;
					double fpr_rec_prev = pop[iter][0].fpr;
					for (i = 0; i < mege_h; i++)success_ri[i] = 0;
					loc_rec = 0;
					for (i = 0; i < mege_h; i++)pop[iter][i].sr = pop[iter][i].tr = 0;
					sr = 0;
					//for (sr = 1; sr <= n_rec_cycle; sr++)
					do
					{
						sr++;
						//if(gen ==1 && (sr-1)%10==0)printf("C%d ", sr);
						BigMixI(pair_take, pair_all);
						for (k = 0; k < pair_all; k++)
						{
							//if (gen == 1)printf("K%d ", k);
							int kk[2];
							for (m = 0; m < 2; m++)kk[m] = pair_d[pair_take[k]][m];
							for (m = 0; m < 2; m++)
							{
								pop[iter][kk[m]].get_copy(&det2[m], nseq, reg_max);
							}
							double fit_parent_max;
							fit_parent_max = Max(det2[0].fit, det2[1].fit);
							int nsn[2], nsi[2], num[2];
							int r_cy;
							if (gen <= 1)r_cy = rand() % 5;
							else r_cy = 4;
							switch (r_cy)
							{
							case(0):
							{
								if (Reco2_Original(&det2[0], &det2[1], nsi, num) == -1)continue;
								else break;
							}
							case(1):
							{
								if (Reco2_Economic(&det2[0], &det2[1], nsi, num) == -1)continue;
								else break;
							}
							case(2):
							{
								if (Reco2_One_dinucleotide_local(&det2[0], &det2[1], olen, reg_max) == -1)continue;
								else break;
							}
							case(3):
							{
								if (Reco2_One_dinucleotide_full(&det2[0], &det2[1]) == -1)continue;
								else break;
							}
							case(4):
							{
								if (Reco2Peak(&det2[0], &det2[1], n_train[iter], xporti) == -1)continue;
								else break;
							}
							}
							//printf("out Rec %d\n",r_cy);
							if (r_cy <= 1)
							{
								for (m = 0; m < 2; m++)
								{
									//	det[m].tot[nsi[m]].print_all();
									nsn[m] = det2[m].order(nsi[m]);
								}
							}
							double dd_r[2] = { 0, 0 };
							int gom[2] = { 0, 0 };
							{
								for (m = 0; m < 2; m++)
								{
									for (int t = 0; t < mege_h; t++)
									{
										if (t == pair_d[k][m])continue;
										gom[m] = GomTown(det2[m], pop[iter][t], n_train[iter], xporti, check_lpd);
										if (gom[m] == -1)
										{
											break;
										}
									}
									if (gom[m] == -1)break;
								}
							}
							if (gom[0] != 1 && gom[1] != 1)
							{
								for (m = 0; m < 2; m++)
								{
									det2[m].fit = EvalMahFIT(&det2[m], n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len);//hoxa_wei, qps, frp, octa_rat,
								}
								double fit_det_max = Max(det2[0].fit, det2[1].fit);
								step_rtry[r_cy]++;
								for (m = 0; m < 2; m++)pop[iter][kk[m]].tr++;
								if (fit_det_max > fit_parent_max)
								{
									for (m = 0; m < 2; m++)pop[iter][kk[m]].sr++;
									int kk_min = Min(kk[0], kk[1]);
									success_ri[kk_min]++;
									success_r1[r_cy]++;
									step_rsuccess[r_cy]++;
									if (fit_det_max > fit_rec_prev)
									{
										loc_rec++;
										rec_jump[kk_min]++;
									}
									success_r_cycle++;
									for (m = 0; m < 2; m++)
									{
										det2[m].get_copy(&pop[iter][kk[m]], nseq, reg_max);
									}
								}
							}
						}
						success_r = 0;
						for (i = 0; i < mege_h; i++)if (success_ri[i] > 0)success_r++;
						if (sr % sr_step_cycle == 0)
						{
							int stepr;
							if (gen <= 1)
							{
								stepr = step_rtry[0] + step_rtry[1] + step_rtry[2] + step_rtry[3];
								if (stepr > 0)ratio_r[0] = (double)(step_rsuccess[0] + step_rsuccess[1] + step_rsuccess[2] + step_rsuccess[3]) / stepr;
								else ratio_r[0] = 0;
							}
							else ratio_r[0] = 0;
							stepr = step_rtry[4];
							if (stepr > 0)ratio_r[1] = (double)(step_rsuccess[4]) / stepr;
							else ratio_r[1] = 0;
							qsort((void*)(&pop[iter][0]), mege_h, sizeof(pop[iter][0]), compare_pop);
							int tryr = 0, sucr = 0, rec_jump1 = 0;
							for (i = 0; i < mege_h; i++)
							{
								if (rec_jump[i] > 0)
								{
									rec_jump1++;
									sucr += pop[iter][i].sr;
									tryr += pop[iter][i].tr;
								}
								else break;
							}
							double sr0_rate;
							if (tryr > 0) sr0_rate = (double)sucr / tryr;
							else sr0_rate = 0;
							printf("Rec %d: %d %d,%d,%d,%d,%d %d OL %f P %f M %f H %g F %f ", sr * pair_all, success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4], success_r_cycle, ratio_r[0], ratio_r[1], pop[iter][0].mah, pop[iter][0].fpr, pop[iter][0].fit);
							fit_rec_prev = pop[iter][elit_rec].fit;
							loc_rec_tot += loc_rec;
							double ratio_per_cycle = pop[iter][0].fit / fit_rec_prev0 - 1;
							printf("L%d RE %f SR rate %f Rec Jump %d RatioSco %f\n", loc_rec_tot, fit_rec_prev, sr0_rate, rec_jump1, ratio_per_cycle);
							for (k = 0; k <= elit_rec; k++)printf("%d M %f H %f F %f\n", k + 1, pop[iter][k].mah, pop[iter][k].fpr, pop[iter][k].fit);
							for (k = 0; k < mege_h; k++)printf("%d ", rec_jump[k]);
							printf("\n");
							tnow = time(NULL);
							printf("%s", ctime(&tnow));
							for (k = 0; k < mege_h; k++)rec_jump[k] = 0;
							if (loc_rec == 0)break;
							if (ratio_r[0] < ratio_thr_r[0] && ratio_r[1] < ratio_thr_r[1])break;
							if (sr0_rate <= sm0_rate)break;
							if (ratio_per_cycle <= ratio_rec_cycle)break;
							if (sr > n_rec_cycle && rec_jump1 <= elit_rec + 1)break;
							for (m = 0; m < NREC; m++)step_rtry[m] = step_rsuccess[m] = 0;
							success_r_cycle = 0;
							loc_rec = 0;
							fit_rec_prev0 = pop[iter][0].fit;
							//for (m = 0; m < mege_h; m++)pop[iter][m].sr = pop[iter][m].tr = 0;
						}
					} while (loc_rec >= 0);
					//tnow = time(NULL);
					//printf("Out rec %s", ctime(&tnow));
					qsort((void*)(&pop[iter][0]), mege_h, sizeof(pop[iter][0]), compare_pop);
					for (k = 1; k < mege_h; k++)
					{
						if (stop_pi[k] == 1 && fit_pre_rec[k] > pop[iter][k].fit)stop_pi[k] = 0;
					}
					double change_level = pop[iter][0].fit / fit_prev - 1;
					double change_level_rec = pop[iter][0].fit / fit_after_mut - 1;
					double change_level_mut = fit_after_mut / fit_prev - 1;
					if (restart <= 1)
					{
						if (restart == 0)
						{
							printf("Restart Half\n");
							tnow = time(NULL);
							printf("%s", ctime(&tnow));
							/*for (i = 0; i < ELIT; i++)
							{
								printf("%d ", i + 1);
								pop[iter][i].print_all(reg_max, nseq);
							}*/
							restart = 1;
							int half = ELIT / 2;
							for (i = half; i < ELIT; i++)
							{
								pop[iter][i - half].get_copy(&pop[iter][i], nseq, reg_max);
								pop[iter][i].init_rand_part_hoxa(n_train[iter], xporti, 20, olen, len_octa, len, octa_prowb, octa_prows);
								EvalMahFIT(&pop[iter][i], n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len);//qps, frp, octa_rat,
								stop_pi[i] = 0;
							}
						}
						else
						{
							double fall = 0.95;
							if (pop[iter][0].mah > fall * mah_rec_prev && pop[iter][0].fpr > fall * fpr_rec_prev)
							{
								restart = 2;
								printf("Restart Full\n");
								check_lpd = 0;
								stop_pi[0] = 0;
								for (i = 1; i < ELIT; i++)
								{
									if (GomTown2(pop[iter][0], pop[iter][i]) == -1)continue;
									pop[iter][0].get_copy(&pop[iter][i], nseq, reg_max);
									pop[iter][i].init_rand_part_hoxa(n_train[iter], xporti, 20, olen, len_octa, len, octa_prowb, octa_prows);
									EvalMahFIT(&pop[iter][i], n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len);//qps, frp, octa_rat,
								}
							}
						}
						qsort((void*)(&pop[iter][0]), ELIT, sizeof(pop[iter][0]), compare_pop);
					}
					else restart++;
					gen++;
					fit_prev = pop[iter][0].fit;
					//if (change_level<GA_EXIT)gen1++;				
					printf("Gen %d Fit %.5f Rat %.5f RatM %.5f RatR %.5f ", gen, pop[iter][0].fit, change_level, change_level_mut, change_level_rec);
					printf("M %d %d,%d,%d R %d ", success_m, success_o, success_l, success_p, success_r1[4]);
					{
						int sumr = 0;
						for (i = 0; i < mege_h; i++)sumr += success_ri[i];
						if (sumr > 0)for (i = 0; i < mege_h; i++)printf("%.2f ", 100 * (double)success_ri[i] / sumr);
					}
					{
						//	printf("One ");
							//time_t tnow;
							//time(&tnow);
						tnow = time(NULL);
						printf("%s", ctime(&tnow));
						//printf("Two ");
					}
					fit_prev = pop[iter][0].fit;
					if (restart > 2 && change_level < GA_EXIT)
					{
						big_exit1 = 1;
						printf("Go out %d iteration\n", iter + 1);
					}
					if (gen == 1)for (i = 0; i < mege_h; i++)stop_li[i] = stop_oi[i] = 1;
					for (i = 0; i < mege_h; i++)stop_pi[i] = 0;
				} while (big_exit1 == 0);
				{
					i = 0;
					for (k = 0; k < nseq; k++)
					{
						if (xport[k] == 0)xportj[i++] = k;
					}
				}
				EvalMahControl(&pop[iter][0], &pop_ext, nseq, nseqb, n_train[iter], n_cntrl[iter], xporti, xportj, fp_rate, cnt_count, n_any_tot, seq_real, seq_back, olen, len, lenb, dav, dcv, qp, prc);
				//for (k = 0; k < n_cntrl[iter]; k++){ fp_rate[cnt_count] = 0.0001; sco_pos[cnt_count] = 0.9; cnt_count++; }
				n_any_tot += (n_train[iter] + nseqb);
				{
					char name[500];
					for (i = 0;; i++)
					{
						if (file_for[i] == '.') { name[i] = '\0'; break; }
						if (file_for[i] == '\n') { name[i] = '\0'; break; }
						if (file_for[i] == '\0') { name[i] = '\0'; break; }
						name[i] = file_for[i];
					}
					char extmat[20];
					char extmat0[] = "_mat";
					strcpy(extmat, extmat0);
					char file_for1[500];
					strcpy(file_for1, path_out);
					strcat(file_for1, file_for);
					pop[iter][0].fprint_allfi_mat(file_for1, extmat, name, olen, pop_ext.c0, pop_ext.buf, iter, size_start, olen_min0);
				}
				big_exit1 = 1;
			}
			qsort(fp_rate, n_cnt_tot, sizeof(double), compare_qq);
			qsort(prc, n_both_sam, sizeof(prc[0]), compare_qbs);
			FILE* outq;
			memset(file_out_cnt, 0, sizeof(file_out_cnt));
			strcpy(file_out_cnt, path_out);
			strcat(file_out_cnt, file_for);
			strcat(file_out_cnt, add_roc);
			if (olen == olen_min0 && size0 == size_start)
			{
				if ((outq = fopen(file_out_cnt, "wt")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			else
			{
				if ((outq = fopen(file_out_cnt, "at")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			if (size0 == size_start)
			{
				double dtp = 1 / (double)CENT;
				double dtp2 = 1 / (double)nseq / 2;
				fprintf(outq, "\t%.5f", dtp2);
				for (n = 1; n < CENT; n++)
				{
					fprintf(outq, "\t%.3f", n * dtp);
				}
				fprintf(outq, "\n");
			}
			fprintf(outq, "%s_%d_%d", file_for, olen, size0);
			for (n = 0; n < CENT; n++)
			{
				fprintf(outq, "\t%.3e", fp_rate[n_decil[n]]);		///qsd											
			}
			fprintf(outq, "\n");
			fclose(outq);
			printf("\nROC\n");
			memset(file_out_cnt, 0, sizeof(file_out_cnt));
			strcpy(file_out_cnt, path_out);
			strcat(file_out_cnt, file_for);
			strcat(file_out_cnt, add_roc1);
			if (olen == olen_min0 && size0 == size_start)
			{
				if ((outq = fopen(file_out_cnt, "wt")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			else
			{
				if ((outq = fopen(file_out_cnt, "at")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			fprintf(outq, "\tROC_%d_%d\n", olen, size0);
			fprintf(outq, "0\t0\n");
			int tproc_pred = 0;
			double fproc_pred = 0;
			double auc_roc = 0;
			int n_cnt_tot1 = n_cnt_tot - 1;
			for (n = 0; n < n_cnt_tot; n++)
			{
				if (fp_rate[n] > fproc_pred && (n == n_cnt_tot1 || fp_rate[n + 1] > fp_rate[n]))
				{
					int tproc_cur = n + 1;
					double fproc_cur = fp_rate[n];
					if (fproc_cur >= fp2 || n == n_cnt_tot1)fproc_cur = fp2;
					double dauc = (tproc_cur + tproc_pred) * (fproc_cur - fproc_pred) / 2 / n_cnt_tot;
					printf("%d\t%d\t%g\t%g\t%g\n", tproc_cur, tproc_pred, fproc_cur, fproc_pred, dauc);
					fprintf(outq, "%g\t%f\n", fproc_cur, (double)tproc_cur / n_cnt_tot);
					auc_roc += dauc;
					if (fproc_cur >= fp2)break;
					tproc_pred = tproc_cur;
					fproc_pred = fproc_cur;
				}
			}
			fclose(outq);
			if (auc_roc > auc_roc_max)
			{
				auc_roc_max = auc_roc;
				size_selected = size0;
				olen_selected = olen;
				for (n = 0; n < n_cnt_tot + 1; n++)fp_rate_best[n] = fp_rate[n];
			}
			if (auc_roc > auc_roc_len)
			{
				auc_roc_len = auc_roc;
				size_len = size0;
			}
			printf("\nTOP 900 scores\n");
			for (n = 0; n < 900; n++)printf("%.12f\t%d\n", prc[n].q, prc[n].n);
			printf("\nTOP 900 FP rates\n");
			for (n = 0; n < 900; n++)printf("%d\t%g\n", n + 1, fp_rate[n]);
			memset(file_out_cnt, 0, sizeof(file_out_cnt));
			strcpy(file_out_cnt, path_out);
			strcat(file_out_cnt, file_for);
			strcat(file_out_cnt, add_prc1);
			if (olen == olen_min0 && size0 == size_start)
			{
				if ((outq = fopen(file_out_cnt, "wt")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			else
			{
				if ((outq = fopen(file_out_cnt, "at")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			fprintf(outq, "\tPRC_%d_%d\n", olen, size0);
			int tpc = 0, fpc = 0;
			printf("\nPRC\n");
			double auc_pr = 0;//prc
			int tpc_pred = tpc;
			double prec_pred;
			if (prc[0].n == 1)prec_pred = 1;
			else prec_pred = 0;
			int prc_count = 0;
			int n_both_sam1 = n_both_sam - 1;
			for (n = 0; n < n_both_sam; n++)
			{
				if (prc[n].n == 0)fpc++;
				else tpc++;
				if (tpc == tpc_pred)continue;
				int n1 = n + 1;
				if (n == n_both_sam1 || ((prc[n].n == 1 && prc[n1].n == 0) && prc[n].q > prc[n1].q))
				{
					prc_count++;
					double prec = (double)tpc / (tpc + fpc);
					if (prc_count == 1)
					{
						prec_pred = prec;
						fprintf(outq, "0\t%f\n", prec);
					}
					double dauc = (prec + prec_pred) * (tpc - tpc_pred) / 2 / n_cnt_tot;
					auc_pr += dauc;
					printf("%g\t%g\t%d\t%d\t%g\n", prec_pred, prec, tpc_pred, tpc, dauc);
					fprintf(outq, "%f\t%f\n", (double)tpc / n_cnt_tot, prec);
					tpc_pred = tpc;
					prec_pred = prec;
				}
			}
			fclose(outq);
			printf("\n");
			if (auc_pr > auc_pr_max)
			{
				auc_pr_max = auc_pr;
				for (n = 0; n < n_both_sam; n++)
				{
					prc_best[n].n = prc[n].n;
					prc_best[n].q = prc[n].q;
				}
				size_selected2 = size0;
				olen_selected2 = olen;
			}
			if (auc_pr > auc_pr_len)
			{
				auc_pr_len = auc_pr;
				size_len2 = size0;
			}
			memset(file_out_cnt, 0, sizeof(file_out_cnt));
			strcpy(file_out_cnt, path_out);
			strcat(file_out_cnt, file_for);
			strcat(file_out_cnt, add_auc);
			if (olen == olen_min0 && size0 == size_start)
			{
				if ((outq = fopen(file_out_cnt, "wt")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			else
			{
				if ((outq = fopen(file_out_cnt, "at")) == NULL)
				{
					printf("Output file can't be opened!\n");
					exit(1);
				}
			}
			fprintf(outq, "%s\t%d\t%d\t%g\t%g\n", file_for, olen, size0, auc_roc, auc_pr);
			fclose(outq);
		}
		memset(file_out_cnt, 0, sizeof(file_out_cnt));
		strcpy(file_out_cnt, path_out);
		strcat(file_out_cnt, file_for);
		strcat(file_out_cnt, "_len");
		strcat(file_out_cnt, add_auc);
		FILE* outq;
		if (olen == olen_min0)
		{
			if ((outq = fopen(file_out_cnt, "wt")) == NULL)
			{
				printf("Output file can't be opened!\n");
				exit(1);
			}
		}
		else
		{
			if ((outq = fopen(file_out_cnt, "at")) == NULL)
			{
				printf("Output file can't be opened!\n");
				exit(1);
			}
		}
		fprintf(outq, "%s\t%d\t%d\t%g\t", file_for, olen, size_len, auc_roc_len);
		fprintf(outq, "%d\t%g\n", size_len2, auc_pr_len);
		fclose(outq);
	}
	{
		FILE* outq;
		memset(file_out_cnt, 0, sizeof(file_out_cnt));
		strcpy(file_out_cnt, path_out);
		strcat(file_out_cnt, file_for);
		strcat(file_out_cnt, "_best");
		strcat(file_out_cnt, argv[5]);
		strcat(file_out_cnt, add_roc);
		if ((outq = fopen(file_out_cnt, "wt")) == NULL)
		{
			printf("Output file can't be opened!\n");
			exit(1);
		}
		{
			double dtp = 1 / (double)CENT;
			double dtp2 = 1 / (double)n_cnt_tot / 2;
			fprintf(outq, "\t%.5f", dtp2);
			for (n = 1; n < CENT; n++)
			{
				fprintf(outq, "\t%.3f", n * dtp);
			}
			fprintf(outq, "\n");
		}
		fprintf(outq, "%s_%d", file_for, size_selected);
		for (n = 0; n < CENT; n++)
		{
			fprintf(outq, "\t%.3e", fp_rate_best[n_decil[n]]);		///qsd											
		}
		fprintf(outq, "\n");
		fclose(outq);
		memset(file_out_cnt, 0, sizeof(file_out_cnt));
		strcpy(file_out_cnt, path_out);
		strcat(file_out_cnt, file_for);
		strcat(file_out_cnt, "_best");
		strcat(file_out_cnt, argv[5]);
		strcat(file_out_cnt, add_prc);
		if ((outq = fopen(file_out_cnt, "wt")) == NULL)
		{
			printf("Output file can't be opened!\n");
			exit(1);
		}
		/*for (n = 0; n < n_both_sam; n++)
		{
			fprintf(outq, "%d\t%.12f\n", prc_best[n].n, prc_best[n].q);
		}*/
		fprintf(outq, "\t%s_%d_%d\n", file_for, olen_selected2, size_selected2);
		double auc_pr1 = 0;//prc
		{
			int tpc = 0, fpc = 0;
			double prec_pred;
			if (prc_best[0].n == 1)prec_pred = 1;
			else prec_pred = 0;
			int tpc_pred = tpc;
			fprintf(outq, "%g\t%g\n", (double)tpc_pred / n_cnt_tot, prec_pred);
			int prc_count = 0;
			int n_both_sam1 = n_both_sam - 1;
			for (n = 0; n < n_both_sam; n++)
			{
				if (prc_best[n].n == 0)fpc++;
				else tpc++;
				if (tpc == tpc_pred)continue;
				int n1 = n + 1;
				if (n == n_both_sam1 || ((prc_best[n].n == 1 && prc_best[n1].n == 0) && prc_best[n].q > prc_best[n1].q))
				{
					double prec = (double)tpc / (tpc + fpc);
					prc_count++;
					if (prc_count == 1)prec_pred = prec;
					double dauc = (prec + prec_pred) * (tpc - tpc_pred) / 2 / n_cnt_tot;
					auc_pr1 += dauc;
					tpc_pred = tpc;
					prec_pred = prec;
					fprintf(outq, "%g\t%g\n", (double)tpc / n_cnt_tot, prec);
				}
			}
			fprintf(outq, "\n");
			fclose(outq);
		}
		memset(file_out_cnt, 0, sizeof(file_out_cnt));
		strcpy(file_out_cnt, path_out);
		strcat(file_out_cnt, file_for);
		strcat(file_out_cnt, "_best");
		strcat(file_out_cnt, argv[5]);
		strcat(file_out_cnt, add_auc);
		if ((outq = fopen(file_out_cnt, "wt")) == NULL)
		{
			printf("Output file can't be opened!\n");
			exit(1);
		}
		//ROC
		int tproc_pred = 0;
		double fproc_pred = 0;
		double auc_roc = 0;
		int n_cnt_tot1 = n_cnt_tot - 1;
		for (n = 0; n < n_cnt_tot; n++)
		{
			if (fp_rate_best[n] > fproc_pred && (n == n_cnt_tot1 || fp_rate_best[n + 1] > fp_rate_best[n]))
			{
				int tproc_cur = n + 1;
				double fproc_cur = fp_rate_best[n];
				if (fproc_cur >= fp2 || n == n_cnt_tot1)fproc_cur = fp2;
				double dauc = (tproc_cur + tproc_pred) * (fproc_cur - fproc_pred) / 2 / n_cnt_tot;
				auc_roc += dauc;
				printf("%d\t%d\t%g\t%g\t%g\n", tproc_cur, tproc_pred, fproc_cur, fproc_pred, dauc);
				if (fproc_cur >= fp2)break;
				tproc_pred = tproc_cur;
				fproc_pred = fproc_cur;
			}
		}
		fprintf(outq, "%s\t%d\t%d\t%g\t", file_for, olen_selected, size_selected, auc_roc);
		fprintf(outq, "%d\t%d\t%g\n", olen_selected2, size_selected2, auc_pr1);
		fclose(outq);
	}
	for (iter = 0; iter < iteration; iter++)
	{
		for (i = 0; i < MEGE; i++)pop[iter][i].mem_out();
		delete[]pop[iter];
	}
	delete[] pop;
	for (i = 0; i < 2; i++)det2[i].mem_out();
	det1.mem_out();
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] seq_back[k][i];
		}
		delete[] seq_back[k];
	}
	delete[] seq_back;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] seq_real[k][i];
		}
		delete[] seq_real[k];
	}
	delete[] seq_real;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dcv[i];
	}
	delete[] dcv;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dav[i];
	}
	delete[] octa_av;
	delete[] octa1;
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_pro1[i];
	}
	delete[] octa_pro1;

	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prow[i];
	}
	delete[] octa_prow;
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prowb[i];
	}
	delete[] octa_prowb;
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prows[i];
	}
	delete[] octa_prows;
	delete[] octa_pro1p;
	delete[] len;
	delete[] lenb;
	delete[] len_octa;
	delete[] thr_octa;
	delete[] xport;
	delete[] xporti;
	delete[] fp_rate;
	delete[] xportj;
	delete[] fp_rate_best;
	delete[] n_train;
	delete[] n_cntrl;
	delete[] octa_rat;
	delete[] qp;
	delete[] prc;
	delete[] prc_best;
	//delete[] hoxa_wei;
	return 0;
}