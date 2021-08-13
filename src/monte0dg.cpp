#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
//#include  <conio.h>
#include  <math.h>
#include  <time.h>
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
//#define LEVEL 7E-7
#define LEVEL 1E-6

struct corr {
	int r1;
	double u[16];
}*cra;
void DelHole(char *str)
{
	char *hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
void Mix(char *a, char *b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMix0(char *d)
{
	int r;
	int len = strlen(d);
	for (r = 0; r < len - 1; r++) Mix(&d[r], &d[1 + r + (rand() % (len - 1 - r))]);
}
int CheckStr(char *file, char *d, int n)
{
	int i, len, ret;
	len = strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		if (strchr("atgc\n", (int)d[i]) != NULL)continue;
		printf("File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
void EvalLen(char *file, int nseq, int *len, double mo[4])
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0, i;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int n = 0;
	int all = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			TransStr(d);
			if (CheckStr(file, d, n) > 0)
			{
				int lenx = strlen(d);
				len[n] = lenx;
				for (i = 0; i < len[n]; i++)
				{
					char c = d[i];
					switch (c) {
					case 'a': {mo[0]++; all++; break; }
					case 'c': {mo[1]++; all++; break; }
					case 'g': {mo[2]++; all++; break; }
					case 't': {mo[3]++; all++; break; }
					default: break;
					}
				}
			}
			else
			{
				printf("Bad sequence N %d\n%s\n", n + 1, d);
				len[n] = 0;				
				//exit(1);
			}			
			n++;
			if (fl == -1)
			{
				for (i = 0; i < 4; i++)mo[i] /= all;
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
			printf("l:%s\nstrlen(l):%d\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%d\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
void EvalSeq(char *file, int &nseq)
{
	char l[SEQLEN];
	int fl = 0;
	FILE  *in;
	nseq = 0;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	while (fgets(l, sizeof(l), in) != NULL)
	{
		if (*l == symbol)nseq++;
	}
	fclose(in);
}
void DelChar(char *str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
void ReadSeq(char *file, int &nseq, int *len, char **seq)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			DelChar(d, ' ');
			DelChar(d, '\t');
			TransStr(d);
			if (CheckStr(file, d, n) > 0)
			{
				d[len[n]] = '\0';
				strcpy(seq[n], d);
				seq[n][len[n]] = '\0';
				n++;
			}
			if (fl == -1)
			{
				fclose(in);
				nseq = n;
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
			printf("l:%s\nstrlen(l):%d\n", l, strlen(l));
			printf("d:%c\nstrlen(d):%d\n", d[0], strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}

int main(int argc, char *argv[])
{
	int i, j, m, k, n;
	char d[SEQLEN], file[100], file2[100];
	double sost[16][2], u[16], sost1[16][2], mo[4];
	int left[2], right[2], rlen[2], *len;
	char **seq;
	char oli[16][3] = {
"aa","ac","ag","at",
"ca","cc","cg","ct",
"ga","gc","gg","gt",
"ta","tc","tg","tt" };

	if (argc != 4)
	{
		printf("%s -1reg -2fileseq  -3fileout", argv[0]);//-shift_min -shift_max
		exit(1);
	}//8 -1 -1 -1 -1 20 1 sf1_93.txtc sf8d
	srand((unsigned)time(NULL));
	int reg1 = atoi(argv[1]);
	int mix;
	int nseq;
	strcpy(file, argv[2]);
	EvalSeq(file, nseq);
	if (nseq <= 100)mix = 20;
	else
	{
		if (nseq > 2000)mix = 1;
		else
		{
			mix = int(2000 / nseq);
		}
	}
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 4; i++)mo[i] = 0;
	EvalLen(file, nseq, len, mo);
	seq = new char*[nseq];
	for (j = 0; j < nseq; j++)
	{
		seq[j] = new char[len[j] + 1];
		if (seq[j] == NULL) { puts("Out of memory..."); exit(1); }
	}
	ReadSeq(file, nseq, len, seq);
	FILE *out2;
	strcpy(file2, argv[3]);
	if ((out2 = fopen(file2, "wt")) == NULL)
	{
		printf("Output file can't be opened!\n");
		exit(1);
	}
	int r1;
	int n_par = 0;
	left[0] = 0;
	for (r1 = 0; r1 < reg1; r1++)
	{
		right[0] = left[0] + r1;
		right[1] = right[0];
		left[1] = left[0];
		int right_most = Max(right[0], right[1]);
		int left_most = Max(-left[0], -left[1]);
		n_par++;
		//printf("NPAR %d\t[%d;%d]\t[%d;%d]\tL1=%d\n",n_par,left[0]+left_most,right[0]+left_most,left[1]+left_most,right[1]+left_most,r1+1);													
	}
	fprintf(out2, "%s\n%d\n%d\n", argv[2], n_par, reg1);
	for (i = 0; i < 4; i++)fprintf(out2, "%f\n", mo[i]);
	fclose(out2);
	if ((cra = new corr[n_par]) == NULL) { puts("Out of memory..."); exit(1); }
	n_par = 0;
	for (r1 = 0; r1 < reg1; r1++)
	{
		right[0] = left[0] + r1;
		right[1] = right[0];
		left[1] = left[0];
		int right_most = Max(right[0], right[1]);
		int left_most = Max(-left[0], -left[1]);
		rlen[1] = -left[1] + right[1] + 1;
		rlen[0] = -left[0] + right[0] + 1;
		//printf("%d\t%d\n",r1+1,n_par);		
		for (j = 0; j < 16; j++)
		{
			for (i = 0; i < 2; i++)sost1[j][i] = 0;
			u[j] = 0;
		}
		double randm = 0;
		double trace[16];
		int cycle;
		double cc[16], cc2[16], cc1[16], test[16];
		for (i = 0; i < 16; i++)trace[i] = cc[i] = cc2[i] = cc1[i] = test[i] = 0;
		for (cycle = 1;; cycle++)
		{
			for (n = 0; n < nseq; n++)
			{
				memset(d, 0, sizeof(d));
				strcpy(d, seq[n]);
				for (m = 0; m < mix; m++)
				{
					BigMix0(d);
					for (int nt = left_most; nt < len[n] - right_most - 1; nt++)
					{
						for (j = 0; j < 16; j++)
						{
							for (i = 0; i < 2; i++)sost[j][i] = 0;
						}
						for (k = 0; k < 2; k++)
						{
							for (i = nt + left[k]; i <= nt + right[k]; i++)
							{
								for (j = 0; j < 16; j++)
								{
									if (strncmp(&d[i], oli[j], 2) == 0)
									{
										sost[j][k]++;
										break;
									}
								}
							}
							for (j = 0; j < 16; j++)sost[j][k] /= rlen[k];
						}
						for (i = 0; i < 16; i++)
						{
							for (k = 0; k < 2; k++)sost1[i][k] += sost[i][k];
							u[i] += sost[i][0] * sost[i][1];
						}
						for (i = 0; i < 16; i++)trace[i] += sost[i][0] * sost[i][1];
						randm++;
					}
				}
			}
			for (i = 0; i < 16; i++)
			{
				cc[i] = trace[i] - sost1[i][0] * sost1[i][1] / randm;
				cc[i] /= randm;
				//cc/=16;
				cc1[i] += cc[i];
				cc2[i] += cc[i] * cc[i];
			}
			//if(cycle%10==0)printf("%d\n",cycle);
			if (cycle % 100 == 0)
			{
				//test=cycle*sqrt((cc2-cc1*cc1/cycle)/(cycle*(cycle-1)))/fabs(cc1);
				int gom = -1;
				for (i = 0; i < 16; i++)
				{
					test[i] = sqrt((cc2[i] - cc1[i] * cc1[i] / cycle) / (cycle*(cycle - 1)));
					if (test[i] > LEVEL) { gom = i; break; }
				}
				if (gom != -1)
				{

					printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bTEST %8d\t%2d\t%.1e", cycle*mix, gom, test[gom]);
				}
				else break;
			}
		}
		for (i = 0; i < 16; i++)
		{
			for (k = 0; k < 2; k++)sost1[i][k] /= randm;
			u[i] /= randm;
		}
		for (i = 0; i < 16; i++)
		{
			u[i] -= sost1[i][0] * sost1[i][1];
		}
		/*if((out=fopen(file,"at"))==NULL)
		{
			printf("Output file can't be opened!\n");
			exit(1);
		}
		fprintf(out1,"[%d;%d] [%d;%d]\t%d\t",left[0]+left_most,right[0]+left_most,left[1]+left_most,right[1]+left_most,r1+1);
		fprintf(out,"[%d;%d] [%d;%d]\tL1=%d\n",left[0]+left_most,right[0]+left_most,left[1]+left_most,right[1]+left_most,r1+1);
		for(i=0;i<16;i++)
		{
			fprintf(out,"%f\t",u[i]);///mo[i1]/mo[j1]/mo[i2]/mo[j2]
		}
		fprintf(out1,"\n");
		fclose(out);
		printf("\t%d [%d;%d] [%d;%d]\tL1=%d\n",n_par,left[0]+left_most,right[0]+left_most,left[1]+left_most,right[1]+left_most,r1+1);
		for(i=0;i<16;i++)
		{
			printf("%f\t%.1e\t",cc1[i]/cycle,test[i]);
			if((i+1)%4==0)printf("\n");
		}
		*/
		if ((out2 = fopen(file2, "at")) == NULL)
		{
			printf("Output file can't be opened!\n");
			exit(1);
		}
		fprintf(out2, "%d\n", r1 + 1);
		for (i = 0; i < 16; i++)
		{
			cra[n_par].u[i] = u[i];
			fprintf(out2, "%f\n", u[i]);
		}
		fclose(out2);
		n_par++;
		//printf("%d\t%d\n",r1+1,n_par);		
	}
	//	fclose(out1);
	return 1;
}
