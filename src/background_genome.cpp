#define _CRT_SECURE_NO_WARNINGS

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include  <time.h>
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define NCHR 50


struct oliq {
	int num;//nomer
	double sco2;// score
	double sco3;// score
	double sco4;// score
};
struct seqm {
	int num;
	int len;//length
	int nat;// % at
	int don;// ready
	double di[16];
	double tri[64];
	double tet[256];	
};
int compare_len(const void *X1, const void *X2)
{
	struct seqm *S1 = (struct seqm *)X1;
	struct seqm *S2 = (struct seqm *)X2;
	if (S1->len - S2->len > 0)return 1;
	if (S1->len - S2->len < 0)return -1;
	return 0;
}
int compare_sco2(const void *X1, const void *X2)
{
	struct oliq *S1 = (struct oliq *)X1;
	struct oliq *S2 = (struct oliq *)X2;
	if (S1->sco2 - S2->sco2 > 0)return 1;
	if (S1->sco2 - S2->sco2 < 0)return -1;
	return 0;
}
int compare_sco3(const void *X1, const void *X2)
{
	struct oliq *S1 = (struct oliq *)X1;
	struct oliq *S2 = (struct oliq *)X2;
	if (S1->sco3 - S2->sco3 > 0)return 1;
	if (S1->sco3 - S2->sco3 < 0)return -1;
	return 0;
}
int compare_sco4(const void *X1, const void *X2)
{
	struct oliq *S1 = (struct oliq *)X1;
	struct oliq *S2 = (struct oliq *)X2;
	if (S1->sco4 - S2->sco4 > 0)return 1;
	if (S1->sco4 - S2->sco4 < 0)return -1;
	return 0;
}
int compare_num(const void *X1, const void *X2)
{
	struct oliq *S1 = (struct oliq *)X1;
	struct oliq *S2 = (struct oliq *)X2;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	return 0;
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
int CheckStr(char *file, char *d, int n, int print)
{
	int i, len, ret;
	len = strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		if (strchr("atgcATGC", (int)d[i]) != NULL)continue;
		if (print == 1)printf("File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
void EvalSeq(char *file, int &nseq, int olen)
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
	nseq = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = strlen(d);
			int check = CheckStr(file, d, n, 1);
			if (lenx >= olen && check == 1)nseq++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelChar(l, '\n');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelChar(l, '\n');
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
		DelChar(l, '\n');
		strcat(d, l);
	}
}
void EvalLen(char *file, int *len, int olen)
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
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = strlen(d);
			int check = CheckStr(file, d, n, 0);
			if (lenx >= olen && check == 1)len[n++] = lenx;
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
			DelChar(l, '\n');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelChar(l, '\n');
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
		DelChar(l, '\n');
		strcat(d, l);
	}
}
void GetSost(char *d, int word, double *sost)
{
	int i, j, k, i_sost, let;
	char letter0[5] = "acgt", letter[5] = "ACGT";
	if (strchr(letter0, (int)d[0]) != NULL)strcpy(letter, letter0);
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	int lens = strlen(d);
	int size = 1;
	for (k = 0; k < word; k++)size *= 4;
	for (i = 0; i < size; i++)sost[i] = 0;
	lens = lens - word + 1;
	for (i = 0; i < lens; i++)
	{
		i_sost = 0;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k < 4; k++)
			{
				if (d[i + j] == letter[k]) { let = k; break; }
			}
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i_sost]++;
	}
	for (i = 0; i < size; i++)sost[i] /= lens;
}
//get Karlin measure
int GetKarl(char *d, double *freq2, double *freq3, double *freq4)
{
	int i;
	double sta1[4], sta2[16], sta3[64], sta4[256];
	double sta2_xnz[16], sta2_xnnw[16], sta2_xynn[16], sta2_xnzn[16], sta2_nyzn[16], sta2_nynw[16], sta2_nnzw[16];
	double sta3_xynw[64], sta3_xnzw[64], sta3_xyzn[64], sta3_nyzw[64];
	int mo_rel[4] = { 3,2,1,0 };
	int di_rel[16] = { 15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0 };
//char letter[5] = "ACGT";
	int fX2[16] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 };
	int fY2[16] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 };

	//compl
	int tri_rel[64] = {63,47,31,15,59,43,27,11,55,39,23,7,51,35,19,3,62,46,30,14,58,42,26,10,54,38,22,6,50,34,18,2,61,45,29,13,57,41,25,9,53,37,21,5,49,33,17,1,60,44,28,12,56,40,24,8,52,36,20,4,48,32,16,0};
	//fXnZ
	int fXnZ[64] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7,8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11,12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15};
	//fXYn
	int fXYn[64] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15};
	//fnYZ
	int fnYZ[64] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	//fX
	int fX3[64] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	//fY
	int fY3[64] = {0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0};
	//fZ
	int fZ3[64] = {0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0};
	
	//compl
	int tet_rel[256]= {255, 191, 127, 63, 239, 175, 111, 47, 223, 159, 95, 31, 207, 143, 79, 15, 251, 187, 123, 59, 235, 171, 107, 43, 219, 155, 91, 27, 203, 139, 75, 11, 247, 183, 119, 55, 231, 167, 103, 38, 215, 151, 87, 23, 199, 135, 71, 7, 243, 179, 115, 48, 227, 163, 99, 35, 211, 147, 83, 19, 195, 131, 67, 3, 254, 190, 126, 62, 238, 174, 110, 46, 222, 158, 94, 30, 206, 142, 69, 14, 250, 186, 122, 58, 234, 170, 106, 42, 218, 154, 78, 26, 202, 138, 74, 10, 246, 182, 118, 54, 230, 166, 86, 38, 214, 150, 86, 22, 198, 134, 70, 6, 242, 178, 93, 50, 226, 162, 98, 34, 210, 146, 82, 18, 194, 130, 66, 2, 253, 189, 125, 61, 237, 173, 109, 45, 221, 157, 93, 29, 205, 107, 77, 13, 249, 185, 121, 57, 233, 169, 105, 41, 217, 113, 89, 25, 201, 137, 73, 9, 245, 181, 117, 53, 229, 118, 101, 37, 213, 149, 85, 21, 197, 133, 69, 5, 241, 122, 113, 49, 225, 161, 97, 33, 209, 145, 81, 17, 193, 129, 65, 1, 252, 188, 124, 60, 236, 172, 108, 44, 220, 156, 92, 28, 129, 140, 76, 12, 248, 184, 120, 56, 232, 168, 104, 40, 132, 152, 88, 24, 200, 136, 72, 8, 244, 180, 116, 52, 134, 164, 100, 36, 212, 148, 84, 20, 196, 132, 68, 4, 135, 176, 112, 48, 224, 160, 96, 32, 208, 144, 80, 16, 192, 128, 64, 0};
	//XnnW
	int fXnnW[256] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15};
	//XY
	int fXYnn[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15};
	//XnZ
	int fXnZn[256] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15};
	//YZ
	int fnYZn[256] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15};
	//YnW
	int fnYnW[256] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15};
	//ZW
	int fnnZW[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	//XYnW
	int fXYnW[256] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 16, 17, 18, 19, 16, 17, 18, 19, 16, 17, 18, 19, 16, 17, 18, 19, 20, 21, 22, 23, 20, 21, 22, 23, 20, 21, 22, 23, 20, 21, 22, 23, 24, 25, 26, 27, 24, 25, 26, 27, 24, 25, 26, 27, 24, 25, 26, 27, 28, 29, 30, 31, 28, 29, 30, 31, 28, 29, 30, 31, 28, 29, 30, 31, 32, 33, 34, 35, 32, 33, 34, 35, 32, 33, 34, 35, 32, 33, 34, 35, 36, 37, 38, 39, 36, 37, 38, 39, 36, 37, 38, 39, 36, 37, 38, 39, 40, 41, 42, 43, 40, 41, 42, 43, 40, 41, 42, 43, 40, 41, 42, 43, 44, 45, 46, 47, 44, 45, 46, 47, 44, 45, 46, 47, 44, 45, 46, 47, 48, 49, 50, 51, 48, 49, 50, 51, 48, 49, 50, 51, 48, 49, 50, 51, 52, 53, 54, 55, 52, 53, 54, 55, 52, 53, 54, 55, 52, 53, 54, 55, 56, 57, 58, 59, 56, 57, 58, 59, 56, 57, 58, 59, 56, 57, 58, 59, 60, 61, 62, 63, 60, 61, 62, 63, 60, 61, 62, 63, 60, 61, 62, 63};
	//XnZW
	int fXnZW[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
	//XYZ
	int fXYZn[256] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33, 34, 34, 34, 34, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 38, 38, 39, 39, 39, 39, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, 62, 62, 63, 63, 63, 63};
	//YZW
	int fnYZW[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
	//fX
	int fX4[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//fY
	int fY4[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//fZ
	int fZ4[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//fW
	int fW4[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	for (i = 0; i < 4; i++)sta1[i] = 0;
	for (i = 0; i < 16; i++)sta2[i] = 0;
	for (i = 0; i < 64; i++)sta3[i] = 0;
	for (i = 0; i < 256; i++)sta4[i] = 0;
	GetSost(d, 1, sta1);
	GetSost(d, 2, sta2);
	GetSost(d, 3, sta3);
	GetSost(d, 4, sta4);
	for (i = 0; i < 2; i++)
	{
		sta1[i] = sta1[mo_rel[i]] = (sta1[i] + sta1[mo_rel[i]]) / 2;
	}
	for (i = 0; i <= 8; i++)
	{
		if (i >= di_rel[i])continue;
		sta2[i] = sta2[di_rel[i]] = (sta2[i] + sta2[di_rel[i]]) / 2;
	}
	for (i = 0; i < 16; i++) freq2[i] = sta2[i] / (sta1[fX2[i]] * sta1[fY2[i]]);
	for (i = 0; i < 64; i++)
	{
		if (i > tri_rel[i])continue;
		sta3[i] = sta3[tri_rel[i]] = (sta3[i] + sta3[tri_rel[i]]) / 2;
	}
	for (i = 0; i < 256; i++)
	{
		if (i >= tet_rel[i])continue;
		sta4[i] = sta4[tet_rel[i]] = (sta4[i] + sta4[tet_rel[i]]) / 2;
	}
	for (i = 0; i < 16; i++) sta2_xnz[i] = sta2_xnnw[i]= sta2_xynn[i] = sta2_xnzn[i] = sta2_nyzn[i] = sta2_nynw[i] = sta2_nnzw[i] = 0;
	for (i = 0; i < 64; i++) sta3_xynw[i] = sta3_xnzw[i] = sta3_xyzn[i] = sta3_nyzw[i] = 0;
	for (i = 0; i < 64; i++)
	{
		sta2_xnz[fXnZ[i]] += sta3[i];
	}
	for (i = 0; i < 256; i++)
	{
		double s4 = sta4[i];
		sta2_xnnw[fXnnW[i]] += s4;
		sta2_xynn[fXYnn[i]] += s4;
		sta2_xnzn[fXnZn[i]] += s4;
		sta2_nyzn[fnYZn[i]] += s4;
		sta2_nynw[fnYnW[i]] += s4;
		sta2_nnzw[fnnZW[i]] += s4;
		sta3_xynw[fXYnW[i]] += s4;
		sta3_xnzw[fXnZW[i]] += s4;
		sta3_nyzw[fnYZW[i]] += s4;
		sta3_xyzn[fXYZn[i]] += s4;
	}
	for (i = 0; i < 64; i++)
	{		
		double z = sta2_xnz[fXnZ[i]] * sta2[fXYn[i]] * sta2[fnYZ[i]];
		if (z > 0)
		{
			freq3[i] = sta3[i] * sta1[fX3[i]] * sta1[fY3[i]] * sta1[fZ3[i]] / z;
		}
		else freq3[i] = 1;
	}
	for (i = 0; i < 256; i++)
	{				
		double z = sta3_xyzn[fXYZn[i]] * sta3_xynw[fXYnW[i]] * sta3_xnzw[fXnZW[i]] * sta3_nyzw[fnYZW[i]] * sta1[fX4[i]] * sta1[fY4[i]] * sta1[fZ4[i]] * sta1[fW4[i]];
		if (z > 0)
		{
			freq4[i] = sta4[i] * sta2_xynn[fXYnn[i]] * sta2_xnzn[fXnZn[i]] * sta2_xnnw[fXnnW[i]] * sta2_nyzn[fnYZn[i]] * sta2_nynw[fnYnW[i]] * sta2_nnzw[fnnZW[i]];
			freq4[i] /= z;
		}
		else freq4[i] = 1;
	}
	return(1);
}
char *TransStrBack(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c >= 97) d[i] = char(c - 32);
		//else break;
	}
	return(d);
}
int ComplStr(char *d)
{
	int i, len;
	len = strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'A': {d[i] = 'T'; break; }
		case 'T': {d[i] = 'A'; break; }
		case 'C': {d[i] = 'G'; break; }
		case 'G': {d[i] = 'C'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
void ReadSeq(char *file, int nseq, int *len, char ***peak_real, int olen)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, j;	
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = strlen(d[0]);
			int check = CheckStr(file, d[0], n, 0);
			nn++;
			if (lenx >= olen && check == 1)
			{
				TransStrBack(d[0]);									
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				for (j = 0; j < 2; j++)
				{
					strcpy(peak_real[j][n], d[j]);
					peak_real[j][n][len[n]] = '\0';
				}
				n++;
			}
			else
			{
				if (lenx < olen)printf("Short peak %d (Len %d) ignored\n", n + 1, lenx);
				if (check == -1)printf("Unusual symbol, peak %d ignored\n%s\n", n + 1, d[0]);
			}
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelChar(l, '\n');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d[0], 0, sizeof(d[0]));
			DelChar(l, '\n');
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
		DelChar(l,'\n');
		strcat(d[0], l);
	}
}
int main(int argc, char *argv[])
{
	int i, j, k;
	char d[SEQLEN], d1[SEQLEN], filesta[10], fileend[10], genome[10];
	char filei[500], fileo1[500], fileo2[500], fileo3[500], fileo4[500], filechr[NCHR][500], path_fasta[500];
	FILE *out, *in_seq[NCHR], *out2, *out3, *out4;
	if (argc != 14)
	{
		puts("Sintax: 1 path_genome 2file in_fa, 3file out_fa 456files out_fasta_DI_TRI_TETRA 7int height 8double mono prec 9-11double fract_di_tri_tetra 12int back_iter 13 char genome (at10 hg38 mm10)");
		exit(1);
	}
	strcpy(filesta, "chr");
	strcpy(fileend, ".plain");
	char name_chr[NCHR][10];

	//human
	char name_chr_hg[24][3] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	int n_chr_hg=24;
	//hg19
//	int sizelo[24]={249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566};
	//hg38
	int sizelo_hg38[24]={248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415};
	//mouse
	char name_chr_mm[21][3] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y"};
	int n_chr_mm = 21;
	//int n_chr=21;
	//mm9	
	//int sizelo[21]={197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555};	
	//mouse mm10
	int sizelo_mm10[21] = {195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698};
	//arabidopsis
	char name_chr_at[5][3] = { "1","2","3","4","5" };
	int sizelo_at10[5] = { 30427671, 19698289, 23459830, 18585056, 26975502 };
	int n_chr_at = 5;
	//Caenorhabditis elegans WBcel235
	char name_chr_ce[6][5] = { "I","II","III","IV","V","X" };
	int sizelo_ce235[6] = { 15072434, 15279421, 13783801, 17493829, 20924180, 17718942 };
	int n_chr_ce = 6;
	//Drosophila melanogaster	
	char name_chr_dm[7][3] = {"2R","2L","3R","3L","X","Y","4"};	
	//dm5
	//int sizelo[7]={21146708,23011544,27905053,24543557,22422827,1351857};
	//int n_chr=6;
	//dm6
	int sizelo_dm6[7]={25286936,23513712,32079331,28110227,23542271,3667352,1348131};	
	int n_chr_dm=7;
	//yeast Saccharomyces cerevisiae R64-1-1
	int sizelo_sc64[16] = {230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066};
	int n_chr_sc=16;
	char name_chr_sc[16][5] = {"I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"};
	//zea mays
//	int sizelo[10]={301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204};
//	char name_chr[10][3] = {"1","2","3","4","5","6","7","8","9","10"};
	//int n_chr=10;
	int sizelo1[NCHR], sizelo2[NCHR];//lengths in bp, Mb
	int n_chr;
	strcpy(path_fasta, argv[1]);
	strcpy(filei, argv[2]);//in_file
	strcpy(fileo1, argv[3]);//out_file
	strcpy(fileo2, argv[4]);//out_file
	strcpy(fileo3, argv[5]);//out_file
	strcpy(fileo4, argv[6]);//out_file
	int height = atoi(argv[7]);
	double mono_prec= atof(argv[8]);//best by Karlin measure 0.05
	double di_fract = atof(argv[9]);//best by Karlin measure 0.05
	double tri_fract = atof(argv[10]);//best by Karlin measure 0.05
	double tetra_fract = atof(argv[11]);//best by Karlin measure 0.05
	int back_iter = atoi(argv[12]);
	strcpy(genome, argv[13]);
	if (strcmp(genome, "at10") == 0)
	{
		n_chr = n_chr_at;
		for (i = 0; i < n_chr; i++)
		{
			sizelo1[i] = sizelo_at10[i];
			strcpy(name_chr[i], name_chr_at[i]);
		}
	}
	else
	{
		if (strcmp(genome, "hg38") == 0)
		{
			n_chr = n_chr_hg;
			for (i = 0; i < n_chr; i++)
			{
				sizelo1[i] = sizelo_hg38[i];
				strcpy(name_chr[i], name_chr_hg[i]);
			}
		}
		else
		{
			if (strcmp(genome, "mm10") == 0)
			{
				n_chr = n_chr_mm;
				for (i = 0; i < n_chr; i++)
				{
					sizelo1[i] = sizelo_mm10[i];
					strcpy(name_chr[i], name_chr_mm[i]);
				}
			}
			else
			{
				if (strcmp(genome, "dm6") == 0)
				{
					n_chr = n_chr_dm;
					for (i = 0; i < n_chr; i++)
					{
						sizelo1[i] = sizelo_dm6[i];
						strcpy(name_chr[i], name_chr_dm[i]);
					}
				}
				else
				{
					if (strcmp(genome, "ce235") == 0)
					{
						n_chr = n_chr_ce;
						for (i = 0; i < n_chr; i++)
						{
							sizelo1[i] = sizelo_ce235[i];
							strcpy(name_chr[i], name_chr_ce[i]);
						}
					}
					else
					{
						if (strcmp(genome, "sc64") == 0)
						{
							n_chr = n_chr_sc;
							for (i = 0; i < n_chr; i++)
							{
								sizelo1[i] = sizelo_sc64[i];
								strcpy(name_chr[i], name_chr_sc[i]);
							}
						}
					}
				}
			}
		}
	}
	int mil = 1000000;
	if (strcmp(genome, "sc64") == 0)mil = 100000;
	for (i = 0; i < n_chr; i++)sizelo2[i] = sizelo1[i] / mil;
	double stop_thr[4] = { 0.99, 0.95, 0.75, 0.67 };// fraction of peaks 100% covered with height background sequences
	int win_gomol = 50;//to compare homology and min paek length
	srand((unsigned)time(NULL));
	for (i = 0; i < n_chr; i++)
	{
		strcpy(filechr[i],path_fasta);
		strcat(filechr[i], filesta);		
		strcat(filechr[i], name_chr[i]);
		strcat(filechr[i], fileend);
		if ((in_seq[i] = fopen(filechr[i], "rt")) == NULL)
		{
			printf("Input file %s can't be opened!\n", filechr[i]);
		}
	}
	int tot_len = 0;
	for (i = 0; i < n_chr; i++)tot_len += sizelo2[i];
	int *len, nseq = 0, olen = win_gomol;
	EvalSeq(filei, nseq, olen);
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	//	peak_wei = new int[nseq];
	//	if (peak_wei == NULL){ puts("Out of memory..."); exit(1); }
	int dnseq = 0;
	EvalLen(filei, len, olen);
	int **hei;
	hei = new int*[nseq];
	if (hei == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < nseq; i++)
	{
		hei[i] = new int[4];
		if (hei[i] == NULL) { puts("Out of memory..."); exit(1); }
	}
	char ***peak_real;
	peak_real = new char**[2];
	if (peak_real == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		peak_real[i] = new char*[nseq];
		for (j = 0; j < nseq; j++)
		{
			peak_real[i][j] = new char[len[j] + 1];
			if (peak_real[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	ReadSeq(filei, nseq, len, peak_real, win_gomol);
	int len_max, len_min;
	seqm *sort;
	sort = new seqm[nseq];
	if (sort == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < nseq; i++)
	{
		double fr2[16], fr3[64], fr4[256];
		int mono = 0;
		for (j = 0; j < len[i]; j++)
		{
			if ((peak_real[0][i][j] == 'A' || peak_real[0][i][j] == 'T') || (peak_real[0][i][j] == 'a' || peak_real[0][i][j] == 't'))mono++;
		}
		sort[i].num = i;
		sort[i].don = 0;
		sort[i].len = len[i];
		sort[i].nat = mono;
		GetKarl(peak_real[0][i], fr2, fr3, fr4);
		for (j = 0; j < 16; j++)sort[i].di[j] = fr2[j];
		for (j = 0; j < 64; j++)sort[i].tri[j] = fr3[j];
		for (j = 0; j < 256; j++)sort[i].tet[j] = fr4[j];
	}
	for (i = 0; i < nseq; i++)for (j = 0; j < 4; j++)hei[i][j] = 0;
	qsort((void*)(&sort[0]), nseq, sizeof(sort[0]), compare_len);
	len_max = sort[nseq-1].len;
	len_min = sort[0].len;		
	int pr_tot = 0;
	int fl = 0;
	int trys = nseq * back_iter/5;
	int iter = 0;	
	int height0 = 10;
	int nseqb = nseq * height0;
	oliq *sele;
	sele = new oliq[nseqb];
	if (sele == NULL) { puts("Out of memory..."); exit(1); }
	int nseqb1 = nseqb - 1;
	int good = 0;
	int gomol = 0;
	while (iter < trys)
	{
		int rr = rand();
		int chr_z = 1, sum = 0, ra = rr % tot_len;
		for (i = 0; i < n_chr; i++)
		{
			sum += sizelo2[i];
			if (ra < sum)
			{
				chr_z = i;
				break;
			}
		}
		int z_len = sizelo2[chr_z];
		rr = rand();
		int rb = rr % z_len;
		rb *= mil;
		rr = rand();
		int rb1 = rr % 1000;
		rb1 *= 1000;
		rb += rb1;
		rr = rand();
		int rb2 = rr % 1000;
		rb += rb2;
		fseek(in_seq[chr_z], (long)(rb), SEEK_SET);
		int check = -1;				
		iter++;
		while (fgets(d, len_max, in_seq[chr_z]) != NULL)
		{
			check = CheckStr(filechr[chr_z], d, 0, 0);
			if (check == 1)break;
		}
		TransStrBack(d);
		if (check == 1)
		{
			int gom = 0;
			for (i = 0; i < nseq; i++)
			{
				for (k = 0; k < 2; k++)
				{
					for (j = 0; j < len[i] - win_gomol + 1; j++)
					{
						if (strncmp(d, &peak_real[k][i][j], win_gomol) == 0)
						{
							gom = 1;
							break;
						}
					}
					if (gom == 1)break;
				}
				if (gom == 1)break;
			}
			if (gom == 1)
			{
				gomol++;					
			}
			else
			{
				int cat = 0;
				for (i = 0; i < len_min; i++)if ((d[i] == 'A' || d[i] == 'T') || (d[i] == 'a' || d[i] == 't')) cat++;
				int len_cur = len_min;
				int done = 0;
				for (i = 0; i < nseq; i++)
				{
					if (sort[i].don >= height0)continue;					
					for (j = len_cur; j < sort[i].len; j++)
					{
						if ((d[j] == 'A' || d[j] == 'T') || (d[j] == 'a' || d[j] == 't'))cat++;
					}
					double mono = fabs((double)(cat - sort[i].nat)) / sort[i].len;
					if (mono < mono_prec)
					{
						strncpy(d1, d, sort[i].len);
						d1[sort[i].len] = '\0';
						sort[i].don++;
						if (sort[i].don == height0)good++;
						double fr2[16], fr3[64], fr4[256];
						GetKarl(d1, fr2, fr3, fr4);
						sele[pr_tot].num = pr_tot;
						sele[pr_tot].sco4 = sele[pr_tot].sco3 = sele[pr_tot].sco2 = 0;
						for (j = 0; j < 16; j++)sele[pr_tot].sco2 += fabs(fr2[j] - sort[i].di[j]);
						sele[pr_tot].sco2 /= 16;
						for (j = 0; j < 64; j++)sele[pr_tot].sco3 += fabs(fr3[j] - sort[i].tri[j]);
						sele[pr_tot].sco3 /= 64;
						for (j = 0; j < 256; j++)sele[pr_tot].sco4 += fabs(fr4[j] - sort[i].tet[j]);
						sele[pr_tot].sco4 /= 256;
						//fprintf(out, ">chr%s_%d_%d #%d FrAT %.2f\n", name_chr[chr_z], rb, rb + sort[i].len,sort[i].don,100*(double)sort[i].nat/sort[i].len); 
				//		fprintf(out, ">peak%d_%d_Di_%f_Tri_%f_Tetra_%f\n", pr_tot+1, sort[i].don, sele[pr_tot].sco2, sele[pr_tot].sco3, sele[pr_tot].sco4);
					//	for (j = 0; j < sort[i].len; j++)fprintf(out, "%c", d[j]);
					//	fprintf(out, "\n");
						pr_tot++;
						done = 1;							
						break;
					}
					len_cur = sort[i].len;
				}
			}
		}								
		if (iter % 10000 == 0)
		{
			double rat = double(good) / nseq;			
			printf("Iterations %5d\t Nseq_Background %5d\tFraction_Done %5f\tHomol %d\n", iter, pr_tot, rat,gomol);
			if (rat > stop_thr[0])break;
		}
	//	if (pr_tot == nseqb1)break;		
	}
	qsort((void*)(&sele[0]), pr_tot, sizeof(sele[0]), compare_sco2);
	double thr_karl2, thr_karl3, thr_karl4;
	{
		int num_karl;
		if (di_fract == 1)num_karl = pr_tot - 1;
		else num_karl = (int(pr_tot*di_fract));
		thr_karl2 = sele[num_karl].sco2;
	}
	qsort((void*)(&sele[0]), pr_tot, sizeof(sele[0]), compare_sco3);
	{
		int num_karl;
		if (tri_fract == 1)num_karl = pr_tot - 1; 
		else num_karl = (int(pr_tot*tri_fract));
		thr_karl3 = sele[num_karl].sco3;
	}
	qsort((void*)(&sele[0]), pr_tot, sizeof(sele[0]), compare_sco4);
	{
		int num_karl;
		if (tetra_fract == 1)num_karl = pr_tot - 1;
		else num_karl = (int(pr_tot*tetra_fract));
		thr_karl4 = sele[num_karl].sco4;
	}
	printf("Score2 %f\tScore3 %f\tScore4 %f\n", thr_karl2, thr_karl3, thr_karl4);	
	if ((out = fopen(fileo1, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo1);
		exit(1);
	}
	if ((out2 = fopen(fileo2, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo2);
		exit(1);
	}
	if ((out3 = fopen(fileo3, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo3);
		exit(1);
	}
	if ((out4 = fopen(fileo4, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo4);
		exit(1);
	}
	int heis[4] = { 0,0,0,0 };
	int size[4] = { 0,0,0,0 };
	iter = pr_tot = 0;
	trys = nseq * back_iter;
	gomol = 0;
	int stop[4];
	{		
		for (i = 0; i < 4; i++)stop[i] = (int)(stop_thr[i] * nseq);
	}	
	while (iter < trys && heis[3] < nseq)
	{
		int rr = rand();
		int chr_z = 1, sum = 0, ra = rr % tot_len;
		for (i = 0; i < n_chr; i++)
		{
			sum += sizelo2[i];
			if (ra < sum)
			{
				chr_z = i;
				break;
			}
		}
		int z_len = sizelo2[chr_z];
		rr = rand();
		int rb = rr % z_len;
		rb *= 1000000;
		rr = rand();
		int rb1 = rr % 1000;
		rb1 *= 1000;
		rb += rb1;
		rr = rand();
		int rb2 = rr % 1000;
		rb += rb2;
		fseek(in_seq[chr_z], (long)(rb), SEEK_SET);
		int check = -1;		
		iter++;
		while (fgets(d, len_max, in_seq[chr_z]) != NULL)
		{
			check = CheckStr(filechr[chr_z], d, 0, 0);
			if (check == 1)break;
		}
		TransStrBack(d);
		if (check == 1)
		{
			int gom = 0;
			for (i = 0; i < nseq; i++)
			{
				for (k = 0; k < 2; k++)
				{
					for (j = 0; j < len[i] - win_gomol + 1; j++)
					{
						if (strncmp(d, &peak_real[k][i][j], win_gomol) == 0)
						{
							gom = 1;
							break;
						}
					}
					if (gom == 1)break;
				}
				if (gom == 1)break;
			}
			if (gom == 1)
			{
				gomol++;					
			}
			else
			{
				int cat = 0;
				for (i = 0; i < len_min; i++)if ((d[i] == 'A' || d[i] == 'T') || (d[i] == 'a' || d[i] == 't')) cat++;
				int len_cur = len_min;
				int done = 0;
				for (i = 0; i < nseq; i++)
				{
					if (hei[i][3] >= height)continue;					
					for (j = len_cur; j < sort[i].len; j++)
					{
						if ((d[j] == 'A' || d[j] == 'T') || (d[j] == 'a' || d[j] == 't'))cat++;
					}
					double mono = fabs((double)(cat - sort[i].nat)) / sort[i].len;
					if (mono < mono_prec)
					{
						strncpy(d1, d, sort[i].len);
						d1[sort[i].len] = '\0';
						double fr2[16], fr3[64], fr4[256];
						GetKarl(d1, fr2, fr3, fr4);
						double sco[3] = { 0, 0, 0 };
						for (j = 0; j < 16; j++)sco[0] += fabs(fr2[j] - sort[i].di[j]);
						sco[0] /= 16;
						for (j = 0; j < 64; j++)sco[1] += fabs(fr3[j] - sort[i].tri[j]);
						sco[1] /= 64;
						for (j = 0; j < 256; j++)sco[2] += fabs(fr4[j] - sort[i].tet[j]);
						sco[2] /= 256;
						//fprintf(out, ">chr%s_%d_%d #%d FrAT %.2f\n", name_chr[chr_z], rb, rb + sort[i].len,sort[i].don,100*(double)sort[i].nat/sort[i].len); 
				//		DelChar(d, '\n');
					//	DelChar(d, '\r');
						if (hei[i][0] < height)
						{
							hei[i][0]++;
							fprintf(out, ">peak%d_%d_Mo_%f_Di_%f_Tri_%f_Tetra_%f\n", sort[i].num, hei[i][0], mono, sco[0], sco[1], sco[2]);
							fprintf(out, "%s\n", d1);
							if (hei[i][0] == height)heis[0]++;
							done = 1;
							size[0]++;
						}
						if (sco[0] < thr_karl2)
						{
							if (hei[i][1] < height)
							{
								hei[i][1]++;
								fprintf(out2, ">peak%d_%d_Mo_%f_Di_%f_Tri_%f_Tetra_%f\n", sort[i].num, hei[i][1], mono, sco[0], sco[1], sco[2]);
								fprintf(out2, "%s\n", d1);
								if (hei[i][1] == height)heis[1]++;
								done = 1;
								size[1]++;
							}
							if (sco[1] < thr_karl3)
							{
								if (hei[i][2] < height)
								{
									hei[i][2]++;
									fprintf(out3, ">peak%d_%d_Mo_%f_Di_%f_Tri_%f_Tetra_%f\n", sort[i].num, hei[i][2], mono, sco[0], sco[1], sco[2]);
									fprintf(out3, "%s\n", d1);
									if (hei[i][2] == height)heis[2]++;
									done = 1;
									size[2]++;
								}
								if (sco[2] < thr_karl4)
								{
									if (hei[i][3] < height)
									{
										hei[i][3]++;
										fprintf(out4, ">peak%d_%d_Mo_%f_Di_%f_Tri_%f_Tetra_%f\n", sort[i].num, hei[i][3], mono, sco[0], sco[1], sco[2]);
										fprintf(out4, "%s\n", d1);
										if (hei[i][3] == height)heis[3]++;
										done = 1;
										size[3]++;
									}
								}
							}
						}
					}
					len_cur = sort[i].len;
					if (done == 1)
					{
						pr_tot++;	
						break;
					}										
				}					
			}					
		}		
		if (iter % 10000 == 0)
		{					
			printf("Iterations %5d\t Nseq_Background %5d\tFraction_Done %5f %5f %5f %5f\tSizes %d %d %d %d\tHomol %d\n", iter, pr_tot, (double)heis[0] / nseq, (double)heis[1] / nseq, (double)heis[2] / nseq,(double)heis[3]/nseq,size[0],size[1],size[2],size[3],gomol);
			int exiti = 1;
			for (i = 0; i < 4; i++)
			{
				if (heis[i] < stop[i])
				{
					exiti = 0;
					break;
				}
			}
			if (exiti == 1)break;
		}
	}
	//qsort((void*)(&sele[0]), pr_tot, sizeof(sele[0]), compare_num);
	for (i = 0; i < n_chr; i++)fclose(in_seq[i]);
	fclose(out);
	fclose(out2); 
	fclose(out3);
	fclose(out4);
	for (j = 0; j < nseq; j++)delete[] hei[j];
	delete[] hei;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] peak_real[k][i];
		}
		delete[] peak_real[k];
	}
	delete[] peak_real;
	delete[] len;	
	delete[] sort;
	delete[] sele;
	return 0;
}