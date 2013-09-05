// HIFI.cpp : Defines the entry point for the console application.
//
//#include "stdafx.h"

#include "stdlib.h"
#include "malloc.h"
#include "string.h"
#include "time.h"
#include "stdio.h"
#include "errno.h"


#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),(mode)))==NULL
#endif


#define LINECAPACITY 10000
#define RSCAPACITY 40
//#define MAFSTEP 0.01	ori minor allele frequence step
//#define MAFSTEP 0.1
#define IMPUTORTABLEWIDTH 1024
#define MAXCOMPARELEN 1000
#define WINI 16
//WINI should be > WINO for the ini
#define WINO 10
#define HEADOFHAPLOLENGTH  1000

int *array_union(int *a, int na, int *b, int nv, int *nc);
int cmp_rise(const void *a, const void *b);
bool hogoxn_cmp(char **a, char **b, int *c, int nc);
void hogoxn_cpy(char **a, char **b, int *c, int nc);
int match_hr(char **hr, int n_hr, char **hxn, int l_hxn, char **hr_refine, char**hr_refine_f, int *hxnmat_0, int *hxnmat_1,int *matind_0, int *matind_1);
int match_hr0(char **hr, int n_hr, char **hxn, int l_hxn, char **hr_refine, char**hr_refine_f, int *hxnmat_0, int *hxnmat_1,int *matind_0, int *matind_1);
int cmporder(const void *a, const void *b);
void cpytable_hr(char **a, char **b, int j);
void setxntohr(char **a, char *b,int j, int k);
int cmphxnhr(char **a, char *b, int j, int *c);
bool matchtest(char *a_1, char *a_2, char **b, int l);
int getHapFre(char *a, char**b,int len);
void getHxnInf(char** a, int len);
void getQScore(char xn, char dsr, int fre_h0, int fre_h1, int fre_total, int winsize_snp, int winsize_pos, double maf_cycno, int haps_sum,int win_homo, int win_hete,int jump_sum, int edge_dist, int win_xx);

//===============================================================\
//Added on 03/20/2012
//1)
static int imp_hap_0 = 0;
static int imp_hap_1 = 0;
static int imp_hap_total = 0;
//2)
static int imp_haps_sum = 0;

static int imp_win_homo = 0;
static int imp_win_hete = 0;
static int imp_win_xx = 0;
static int imp_win_nn = 0;

static double	q_score_haplo = 0;
static double	q_score_geno = 0;

//===============================================================


//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char** argv)
{
	
/*	if(argc!=4){
		printf("The parameter number was wrong!\n");
		return 0;
	}
*/

	char haplotype[32];
	char genotype[32];
	char refhap[32];
	double MAFSTEP = 0.1;
		
	if(argc == 1){
		strcpy(haplotype,"haplotype.txt");
		strcpy(genotype,"genotype.txt");
		strcpy(refhap,"refHaplos.txt");
	}
	
	else if(argc == 2){
		strcpy(haplotype,argv[1]);
		strcpy(genotype,"genotype.txt");
		strcpy(refhap,"refHaplos.txt");
	}
	
	else if(argc == 4){
		strcpy(haplotype,argv[1]);
		strcpy(genotype,argv[2]);
		strcpy(refhap,argv[3]);
	}
	
	else if(argc == 5){
		strcpy(haplotype,argv[1]);
		strcpy(genotype,argv[2]);
		strcpy(refhap,argv[3]);
		MAFSTEP = atof(argv[4]);	// alternative strtod
	}
	else {
		printf("The parameter number was wrong!\n");
		return 0;
	}

	char impute[32];
	char qscore[32];

	strcpy(impute,"imputed_");
	strncat(impute,haplotype,32);
	
	strcpy(qscore,"qscore_");
	strncat(qscore,haplotype,32);
		
	//==========================
	//Runtime caculation
	clock_t start, finish;
	double duration;
	start = clock();
	//==========================

	int imp_double_side_match = 0; //impute double
	int imp_single_side_match = 0;
	int imp_random_maf_cnt = 0;
	int imp_random_maf_notcnt = 0;

	int jump_sum = 0;

//	FILE *fpo_l;
//	fpo_l=fopen("match_length.txt","w");

	FILE *fpo;
	int err;
	//errno_t err;
	FILE *fp;
	int snpsum_h = 0;
	int snpsum_g = 0;
	int snpsum_r = 0;
	int snp_no = 0;
	char line[LINECAPACITY];
	char headofhaplo[HEADOFHAPLOLENGTH];

	printf("**********   Starting... please wait   **********\n");
//1> open haplotype file
	if((err = fopen_s(&fp,haplotype,"rt"))!=NULL)
	{
		printf("Can not open haplotype file!");
		return 0;
	}
	
	fgets(line,LINECAPACITY,fp);
	strcpy(headofhaplo, line);
	
	while(feof(fp) == 0)
	{
		fgets(line,LINECAPACITY,fp);
		snpsum_h++;
	}

	//Allocate the memory for rs & position & allele for haplotype
   char **h_rs = (char**)malloc(snpsum_h*sizeof(char*));
   for(int i=0;i<snpsum_h;i++)
		h_rs[i]=(char*)malloc(RSCAPACITY*sizeof(char));

   int *h_pos = (int*)malloc(snpsum_h*sizeof(int));

   char *h_allele = (char*)malloc(snpsum_h*sizeof(char));

	rewind(fp);
	snp_no = 0;
	fgets(line,LINECAPACITY,fp);
	
	while(fscanf(fp,"%s %d %c\n",h_rs[snp_no],&h_pos[snp_no],&h_allele[snp_no])!=EOF)
	{	
		snp_no++;
	}
	
	fclose(fp);

//2> open genotype file
	if((err = fopen_s(&fp,genotype,"rt"))!=NULL)
	{
		printf("Can not open genotype file!");
		return 0;
	}

	fgets(line,LINECAPACITY,fp);

	while(feof(fp) == 0)
	{
		fgets(line,LINECAPACITY,fp);
		snpsum_g++;
	}

	//Allocate the memory for rs & position & alleles for genotype
   char **g_rs=(char**)malloc(snpsum_g*sizeof(char*));
   for(int i=0;i<snpsum_g;i++)
		g_rs[i]=(char*)malloc(RSCAPACITY*sizeof(char));

   int *g_pos = (int*)malloc(snpsum_g*sizeof(int));

   char **g_alleles=(char**)malloc(snpsum_g*sizeof(char*));
   for(int i=0;i<snpsum_g;i++)
	   g_alleles[i]=(char*)malloc(2*sizeof(char));

   rewind(fp);
   snp_no = 0;
   fgets(line,LINECAPACITY,fp);

   	while(fscanf(fp,"%s %d %c%c\n",g_rs[snp_no],&g_pos[snp_no],&g_alleles[snp_no][0],&g_alleles[snp_no][1])!=EOF)
	{	
		snp_no++;
	}
	
	fclose(fp);	
	
//3> open reference haplotypes file
	if((err = fopen_s(&fp,refhap,"rt"))!=NULL)
	{
		printf("Can not open reference haplotypes file!");
		return 0;
	}

	fgets(line,LINECAPACITY,fp);

	while(feof(fp) == 0)
	{
		fgets(line,LINECAPACITY,fp);
		snpsum_r++;
	}

	//Allocate the memory for rs & position & alleles for reference haplotypes
   char **r_rs=(char**)malloc(snpsum_r*sizeof(char*));
   for(int i=0;i<snpsum_r;i++)
		r_rs[i]=(char*)malloc(RSCAPACITY*sizeof(char));

   int *r_pos = (int*)malloc(snpsum_r*sizeof(int));

    rewind(fp);
	fgets(line,LINECAPACITY,fp);

   char t1[20];
   int t2;
   fscanf(fp,"%s %d",t1,&t2);
   fgets(line,LINECAPACITY,fp);

   int j = 0;
   int ref_num = 0;
   while(line[j] != '\0')
   {
		if((line[j]!='\t')&&(line[j]!=' ')&&(line[j]!='\n'))
			ref_num++;	
	   j++;
   }
 
   char **r_alleles=(char**)malloc(snpsum_r*sizeof(char*));
   for(int i=0;i<snpsum_r;i++)
	   r_alleles[i]=(char*)malloc(ref_num*sizeof(char));

   rewind(fp);
   snp_no = 0;
    fgets(line,LINECAPACITY,fp);

   	while(feof(fp)==0)
	{	
		fscanf(fp,"%s %d",r_rs[snp_no],&r_pos[snp_no]);
		fgets(line,LINECAPACITY,fp);

		j = 0;
		ref_num = 0;
		while(line[j] != '\0')
	   {
			if((line[j]!='\t')&&(line[j]!=' ')&&(line[j]!='\n'))
			{
				r_alleles[snp_no][ref_num] = line[j];
				ref_num++;
			}
		   j++;
	   }

		snp_no++;
	}
	
	fclose(fp);	

	//HoGoMerger
	int *hgunion = NULL, n_hgunion=0;
	hgunion = array_union(h_pos, snpsum_h, g_pos, snpsum_g, &n_hgunion);

	//Allocate hgr_pos
	int *hgr_pos=NULL, snpsum_hgr=0;
	hgr_pos = array_union(hgunion, n_hgunion, r_pos, snpsum_r, &snpsum_hgr);
	free(hgunion);

	//Allocate hgr_rs
	char **hgr_rs=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		hgr_rs[i]=(char*)malloc(RSCAPACITY*sizeof(char));

	//Allocate h_hgrallele
	char *h_hgrallele=(char*)malloc(snpsum_hgr*sizeof(char));

	//Allocate g_hgralleles
	char **g_hgralleles=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		g_hgralleles[i]=(char*)malloc(2*sizeof(char));

	//Allocate r_hgralleles
	char **r_hgralleles=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		r_hgralleles[i]=(char*)malloc(ref_num*sizeof(char));

	//Allocate h_index
	 int *h_index = (int*)malloc(snpsum_h*sizeof(int));

	//Merger
	//1> merger rs, h_hgrallele, g_hgralleles, r_hgralleles
	int m=0, k=0, l=0;
	bool rscopy=false;

	for(int i=0;i<snpsum_hgr;i++)
	{
		rscopy = false;
		if(hgr_pos[i]==h_pos[m])
		{
			if(rscopy==false)
			{
				strcpy(hgr_rs[i],h_rs[m]);
				rscopy=true;
			}
			h_hgrallele[i] = h_allele[m];
			h_index[m]=i;
			m++;
		}
		else
			h_hgrallele[i] = 'N';
		
		if(hgr_pos[i]==g_pos[k])
		{
			if(rscopy==false)
			{
				strcpy(hgr_rs[i],g_rs[k]);
				rscopy=true;
			}
			g_hgralleles[i][0] = g_alleles[k][0];
			g_hgralleles[i][1] = g_alleles[k][1];
			k++;
		}
		else
		{
			g_hgralleles[i][0] = 'N';
			g_hgralleles[i][1] = 'N';
		}

		if(hgr_pos[i]==r_pos[l])
		{
			if(rscopy==false)
			{
				strcpy(hgr_rs[i],r_rs[l]);
				rscopy=true;
			}
			for(int p=0; p<ref_num; p++)
				r_hgralleles[i][p] = r_alleles[l][p];			
			l++;
		}
		else
			for(int p=0; p<ref_num; p++)
				r_hgralleles[i][p] = 'N';	
	}

	//Free h_rs, h_pos, h_allele, g_rs, g_pos, g_alleles, r_rs, r_pos, r_alleles
	//free h_rs
	for(int i=0;i<snpsum_h;i++)
		free(h_rs[i]);
	free(h_rs);
	//free h_pos
	free(h_pos);
	//free h_allele
	free(h_allele);
	//free g_rs
	for(int i=0;i<snpsum_g;i++)
		free(g_rs[i]);
	free(g_rs);
	//free g_pos
	free(g_pos);
	//free g_alleles
	for(int i=0;i<snpsum_g;i++)
		free(g_alleles[i]);
	free(g_alleles);
	//free r_rs
	for(int i=0;i<snpsum_r;i++)
		free(r_rs[i]);
	free(r_rs);
	//free r_pos
	free(r_pos);
	//free r_alleles
	for(int i=0;i<snpsum_r;i++)
		free(r_alleles[i]);
	free(r_alleles);

	//Allocate discrepancy
	char *discrepancy=(char*)malloc(snpsum_hgr*sizeof(char));
	//Allocate array for merge hogo, the known haplotype is hogo_xn[n][0]
	char **hogo_xn=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		hogo_xn[i]=(char*)malloc(2*sizeof(char));
	//Allocate hogo_xn_ori
	char **hogo_xn_ori=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		hogo_xn_ori[i]=(char*)malloc(2*sizeof(char));

	//ini discrepancy
	for(int t=0;t<snpsum_hgr;t++)
		discrepancy[t]='-';

	for(int n=0;n<snpsum_hgr;n++)
	{
		if(((g_hgralleles[n][0]=='N')&&(g_hgralleles[n][1]=='N'))&&(h_hgrallele[n]=='N'))
		{
			hogo_xn[n][0]='N';
			hogo_xn[n][1]='N';
		}
		else if(((g_hgralleles[n][0]!='N')&&(g_hgralleles[n][1]!='N'))&&(h_hgrallele[n]=='N'))
		{
			if(g_hgralleles[n][0]==g_hgralleles[n][1])
			{
				hogo_xn[n][0]=g_hgralleles[n][0];
				hogo_xn[n][1]=g_hgralleles[n][0];
			}
			else
			{
				hogo_xn[n][0]='X';
				hogo_xn[n][1]='X';
			}
		}
		else if(((g_hgralleles[n][0]=='N')&&(g_hgralleles[n][1]=='N'))&&(h_hgrallele[n]!='N'))
		{			
			hogo_xn[n][0]=h_hgrallele[n];
			hogo_xn[n][1]='N';
		}
		else if(((g_hgralleles[n][0]!='N')&&(g_hgralleles[n][1]!='N'))&&(h_hgrallele[n]!='N'))
		{
			if(g_hgralleles[n][0]==g_hgralleles[n][1])
			{
				if(h_hgrallele[n]==g_hgralleles[n][0])
				{
					hogo_xn[n][0]=g_hgralleles[n][0];
					hogo_xn[n][1]=g_hgralleles[n][0];
				}
				else
				{
					hogo_xn[n][0]=g_hgralleles[n][0];
					hogo_xn[n][1]=g_hgralleles[n][0];
					discrepancy[n]='*';
				}
			}
			else
			{
				if((h_hgrallele[n]==g_hgralleles[n][0])||(h_hgrallele[n]==g_hgralleles[n][1]))
				{
					if(h_hgrallele[n]==g_hgralleles[n][0])
					{
						hogo_xn[n][0]=g_hgralleles[n][0];
						hogo_xn[n][1]=g_hgralleles[n][1];
					}
					else
					{
						hogo_xn[n][0]=g_hgralleles[n][1];
						hogo_xn[n][1]=g_hgralleles[n][0];
					}
				}
				else
				{
					hogo_xn[n][0]='X';
					hogo_xn[n][1]='X';
					discrepancy[n]='*';
				}
			}
		}
		else
			printf("Error type 1");
	}

	//Free h_hgrallele
	//g_hgralleles free for later
	free(h_hgrallele);

	//Copy hogo_xn to hogo_xn_ori
	for(int i=0;i<snpsum_hgr;i++)
	{
		hogo_xn_ori[i][0]=hogo_xn[i][0];
		hogo_xn_ori[i][1]=hogo_xn[i][1];
	}

	//Caculate maf and main allele, minor allele
	//Allocate main_allele in main_allele[i][0]
	char **main_allele=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		main_allele[i]=(char*)malloc(2*sizeof(char));
	double *maf = (double*)malloc(snpsum_hgr*sizeof(double));

	char *ralleles_temp=(char*)malloc(ref_num*sizeof(char));

	for(int i=0;i<snpsum_hgr;i++)
	{
		strcpy(ralleles_temp,r_hgralleles[i]);
		qsort(ralleles_temp,ref_num,sizeof(char),cmp_rise);

		if(ralleles_temp[0]==ralleles_temp[ref_num-1])
		{
			main_allele[i][0]=ralleles_temp[0];
			maf[i]=0;
		}
		else
		{
			for(int j=0;j<ref_num-1;j++)
			{
				if(ralleles_temp[j]!=ralleles_temp[j+1])
				{
					double maf_temp = (double)(j+1)/(double)ref_num;
					if(maf_temp > 0.5)
					{
						main_allele[i][0]=ralleles_temp[0];
						main_allele[i][1]=ralleles_temp[ref_num-1];
						maf[i]=1.0-maf_temp;
					}
					else
					{
						main_allele[i][0]=ralleles_temp[ref_num-1];
						main_allele[i][1]=ralleles_temp[0];
						maf[i]=maf_temp;
					}
				}
			}
		}
	}
	free(ralleles_temp);

	//Allocate hogo_xn_pre
	char **hogo_xn_pre=(char**)malloc(snpsum_hgr*sizeof(char*));
	for(int i=0;i<snpsum_hgr;i++)
		hogo_xn_pre[i]=(char*)malloc(2*sizeof(char));

	//Allocate imputor_table
	char **imputor_table=(char**)malloc(ref_num*sizeof(char*));
	for(int i=0;i<ref_num;i++)
		imputor_table[i]=(char*)malloc(IMPUTORTABLEWIDTH*sizeof(char));
	//Allocate imputor_refinetable
	char **imputor_refinetable=(char**)malloc(ref_num*sizeof(char*));
	for(int i=0;i<ref_num;i++)
		imputor_refinetable[i]=(char*)malloc(IMPUTORTABLEWIDTH*sizeof(char));
	//Allocate imputor_refinetable_f
	char **imputor_refinetable_f=(char**)malloc(ref_num*sizeof(char*));
	for(int i=0;i<ref_num;i++)
		imputor_refinetable_f[i]=(char*)malloc(IMPUTORTABLEWIDTH*sizeof(char));
	//Allocate hogoxn_table
	char **hogoxn_table=(char**)malloc(2*sizeof(char*));
	for(int i=0;i<2;i++)
		hogoxn_table[i]=(char*)malloc(IMPUTORTABLEWIDTH*sizeof(char));
	//Allocate badimpsnp_index
	int *badimpsnp_index=(int*)malloc(snpsum_hgr*sizeof(int));
	//Allocate imputeisxn
	char *imputeisxn=(char*)malloc(snpsum_hgr*sizeof(char));
	//Allocate indexofsnp_act
	int *indexofsnp_act=(int*)malloc(IMPUTORTABLEWIDTH*sizeof(int));
	//Allocate hxnmatch_0
	int *hxnmatch_0 = (int*)malloc(ref_num*sizeof(int));
	//Allocate hxnmatch_1
	int *hxnmatch_1 = (int*)malloc(ref_num*sizeof(int));
	//==========================================================================
	//Added on 03/20/2012
	char *imputeisdsr=(char*)malloc(snpsum_hgr*sizeof(char));
	for(int t=0;t<snpsum_hgr;t++)
		imputeisdsr[t]='-';
	//1)
	int *imputed_hap_0=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_hap_0[t]=0;
	int *imputed_hap_1=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_hap_1[t]=0;
	int *imputed_hap_total=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_hap_total[t]=0;
	//2)
	int *imputed_winsize_snp=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_winsize_snp[t]=0;
	//3)
	int *imputed_winsize_pos=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_winsize_pos[t]=0;
	//4)
	double *imputed_maf_cycno=(double*)malloc(snpsum_hgr*sizeof(double));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_maf_cycno[t]=0;
	//5)
	double *imputed_maf_ref=(double*)malloc(snpsum_hgr*sizeof(double));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_maf_ref[t]=0;
	//6)
	int *imputed_haps_sum=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_haps_sum[t]=0;
	//7)
	int *imputed_win_homo=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_win_homo[t]=0;
	int *imputed_win_hete=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_win_hete[t]=0;
	//8)
	int *imputed_jump_sum=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_jump_sum[t]=0;
	//9)
	int *imputed_edge_dist=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_edge_dist[t]=0;
	//10)
	int *imputed_win_xx=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_win_xx[t]=0;
	int *imputed_win_nn=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int t=0;t<snpsum_hgr;t++)
		imputed_win_nn[t]=0;

	//==========================================================================
	//ini imputeisxn
	for(int t=0;t<snpsum_hgr;t++)
		imputeisxn[t]='-';
	bool cmphogoxn_issame=false;
	//Imputing target haplotype (Hxr)
	for(double impute_maf=0.5-MAFSTEP;impute_maf>=0-MAFSTEP;impute_maf=impute_maf-MAFSTEP)
	{
		//Count snp that maf >=impute_maf
		int snpsum_impmaf = 0;
		for(int i=0;i<snpsum_hgr;i++)
		{
			if(maf[i]>=impute_maf)
				snpsum_impmaf++;
		}
		//Allocate imputemaf snps
		int *snpindex_impmaf = (int*)malloc(snpsum_impmaf*sizeof(int));
		int p=0;
		for(int i=0;i<snpsum_hgr;i++)
		{
			if(maf[i]>=impute_maf)
			{
				snpindex_impmaf[p]=i;
				p++;
			}
		}
		//Allocate impute snps base on maf and known snps in haplotype
		int *snpimp=NULL, n_snpimp=0;
		snpimp = array_union(h_index,snpsum_h,snpindex_impmaf,snpsum_impmaf,&n_snpimp);
		//Free snpindex_impmaf
		free(snpindex_impmaf);



		//======================================================================
		//
		//
		//
		//
		//
		//
		//======================================================================
	

		//Impute XX first
		//The NX Presence Monitor				
		do{
			hogoxn_cpy(hogo_xn, hogo_xn_pre, snpimp, n_snpimp);
			//Impute processing
			//NX locator
			for(int s=0;s<n_snpimp;s++)
			{
				if((hogo_xn[snpimp[s]][0]=='X')&&(hogo_xn[snpimp[s]][1]=='X'))
				{	
					 jump_sum = 0;
					bool imputecycle = false;
					bool lessimputemaxlen = false;
					//ini badimpsnp_index
					for(int t=0;t<n_snpimp;t++)
						badimpsnp_index[t]=0;
					//Imputor
					//2> ini win = INIWIN
					int solution_pre1=0,solution_pre2=0;
						int win_half = WINI/2, win_top, win_bottom;
						do{
								win_top=s-win_half;
								win_bottom=s+win_half;
								if(win_top<0)
									win_top=0;
								if(win_bottom>n_snpimp-1)
									win_bottom=n_snpimp-1;

								if(win_bottom-win_top==0)
									break;
								//hogoxn_table
								int u=0;
								for(int i=0;i<2;i++)
								{
									u=0;
									for(int q=0;q<win_bottom-win_top+1;q++)		
									{									
										if(badimpsnp_index[win_top+q]==0)
										{
											hogoxn_table[i][u] = hogo_xn[snpimp[win_top+q]][i];
											indexofsnp_act[u] = win_top+q;
											u++;
										}
									}
									hogoxn_table[i][u] = '\0';
								}

								if(u<MAXCOMPARELEN)
									lessimputemaxlen = true;
								else
									lessimputemaxlen = false;

								//imputor_table
								for(int i=0;i<ref_num;i++)
								{
									u=0;
									for(int q=0;q<win_bottom-win_top+1;q++)		
									{
										if(badimpsnp_index[win_top+q]==0)
										{
											imputor_table[i][u] = r_hgralleles[snpimp[win_top+q]][i];
											u++;
										}
									}
									imputor_table[i][u] = '\0';
								}

								int matindex_0=0,matindex_1=0;
								int solution = match_hr(imputor_table,ref_num, hogoxn_table,u,imputor_refinetable,imputor_refinetable_f,hxnmatch_0,hxnmatch_1,&matindex_0,&matindex_1);						
																
								if((solution_pre1==2)&&(solution_pre2==0))
								{
									jump_sum++;

									if(solution == 0)
									{
										badimpsnp_index[win_top]=1;
										badimpsnp_index[win_bottom]=1;
										solution = 2;
									}
									else if(solution == 2)
									{	
										if(win_top-1>=0)
											badimpsnp_index[win_top-1]=1;
										if(win_bottom+1<n_snpimp)
											badimpsnp_index[win_bottom+1]=1;
									}
								}

								if(solution==1)
								{
									if(u>=WINO)
									{
										for(int i=0;i<u;i++)
										{
											if((hogoxn_table[0][i]=='X')&&(hogoxn_table[1][i]=='X'))
											{
												//====================================
												//added on 03/16/2012
												imp_double_side_match++;
										//		fprintf(fpo_l,"%d\n",u);
												//====================================
												imputeisxn[snpimp[indexofsnp_act[i]]]='x';
												imputeisdsr[snpimp[indexofsnp_act[i]]]='d';
												//====================================
												//Added on 03/20/2012
												imputed_hap_0[snpimp[indexofsnp_act[i]]]=imp_hap_0;
												imputed_hap_1[snpimp[indexofsnp_act[i]]]=imp_hap_1;
												imputed_hap_total[snpimp[indexofsnp_act[i]]]=imp_hap_total;

												imputed_winsize_snp[snpimp[indexofsnp_act[i]]]=u;
												imputed_winsize_pos[snpimp[indexofsnp_act[i]]] = hgr_pos[snpimp[indexofsnp_act[u-1]]] - hgr_pos[snpimp[indexofsnp_act[0]]];
												imputed_maf_cycno[snpimp[indexofsnp_act[i]]] = impute_maf;
												imputed_maf_ref[snpimp[indexofsnp_act[i]]] = impute_maf;
												imputed_haps_sum[snpimp[indexofsnp_act[i]]] = imp_haps_sum;

												imputed_win_homo[snpimp[indexofsnp_act[i]]] = imp_win_homo;
												imputed_win_hete[snpimp[indexofsnp_act[i]]] = imp_win_hete;
												imputed_win_xx[snpimp[indexofsnp_act[i]]] = imp_win_xx;
												imputed_win_nn[snpimp[indexofsnp_act[i]]] = imp_win_nn;

												imputed_jump_sum[snpimp[indexofsnp_act[i]]] = jump_sum;
												imputed_edge_dist[snpimp[indexofsnp_act[i]]] = (u-1-i)>i?i:u-1-i;





												//====================================


						/*fu						printf("Imputing   %s %d   %c%c <-- %c%c\n",
													        hgr_rs[snpimp[indexofsnp_act[i]]],hgr_pos[snpimp[indexofsnp_act[i]]],hogo_xn[snpimp[indexofsnp_act[i]]][0],hogo_xn[snpimp[indexofsnp_act[i]]][1],
															imputor_refinetable[matindex_0][i],  imputor_refinetable[matindex_1][i]);
							*/					
												hogo_xn[snpimp[indexofsnp_act[i]]][0] = imputor_refinetable[matindex_0][i];
												hogo_xn[snpimp[indexofsnp_act[i]]][1] = imputor_refinetable[matindex_1][i];											
											}
										}
									}

									imputecycle = false;
								}
								else if(solution == 0)
								{
									win_half = win_half-1;
									imputecycle = true;
								}
								else if(solution == 2)
								{
									win_half = win_half+2;
									imputecycle = true;
								}
								else
									printf("Error type 3");

							solution_pre1 = solution_pre2;
							solution_pre2=solution;

						}while((imputecycle==true)&&(lessimputemaxlen==true));
				}
			}
			cmphogoxn_issame = hogoxn_cmp(hogo_xn,hogo_xn_pre,snpimp,n_snpimp);
		}while(!cmphogoxn_issame);

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/* fu		fpo=fopen("imputedhaplotype_1.txt","w");	
		fprintf(fpo,"%s", headofhaplo);
		for(int i=0;i<snpsum_hgr;i++)
		{
			if(i==snpsum_hgr-1)
				fprintf(fpo,"%s %d %c %c", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
			else
				fprintf(fpo,"%s %d %c %c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
		}
		fclose(fpo);
*/		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//======================================================================
		//
		//
		//
		//
		//
		//======================================================================

		//Second, impute NN or remain XX
		//The NX Presence Monitor
		do{
			hogoxn_cpy(hogo_xn, hogo_xn_pre, snpimp, n_snpimp);
			//Impute processing
			//NX locator
			for(int s=0;s<n_snpimp;s++)
			{
				if(((hogo_xn[snpimp[s]][0]=='N')||(hogo_xn[snpimp[s]][1]=='N'))||((hogo_xn[snpimp[s]][0]=='X')&&(hogo_xn[snpimp[s]][1]=='X')))
				{		
					//===============================
					jump_sum = 0;
					bool imputecycle = false;
					bool lessimputemaxlen = false;
					//ini badimpsnp_index
					for(int t=0;t<n_snpimp;t++)
						badimpsnp_index[t]=0;
					//Imputor
					//2> ini win = INIWIN
					int solution_pre1=0,solution_pre2=0;
						int win_half = WINI/2, win_top, win_bottom;
						do{
								win_top=s-win_half;
								win_bottom=s+win_half;
								if(win_top<0)
									win_top=0;
								if(win_bottom>n_snpimp-1)
									win_bottom=n_snpimp-1;

								if(win_bottom-win_top==0)
									break;
								//hogoxn_table
								int u=0;
								for(int i=0;i<2;i++)
								{
									u=0;
									for(int q=0;q<win_bottom-win_top+1;q++)		
									{									
										if(badimpsnp_index[win_top+q]==0)
										{
											hogoxn_table[i][u] = hogo_xn[snpimp[win_top+q]][i];
											indexofsnp_act[u] = win_top+q;
											u++;
										}
									}
									hogoxn_table[i][u] = '\0';
								}

								if(u<MAXCOMPARELEN)
									lessimputemaxlen = true;
								else
									lessimputemaxlen = false;

								//imputor_table
								for(int i=0;i<ref_num;i++)
								{
									u=0;
									for(int q=0;q<win_bottom-win_top+1;q++)		
									{
										if(badimpsnp_index[win_top+q]==0)
										{
											imputor_table[i][u] = r_hgralleles[snpimp[win_top+q]][i];
											u++;
										}
									}
									imputor_table[i][u] = '\0';
								}

								int matindex_0=0,matindex_1=0;
								int solution = match_hr(imputor_table,ref_num, hogoxn_table,u,imputor_refinetable,imputor_refinetable_f,hxnmatch_0,hxnmatch_1,&matindex_0,&matindex_1);						
																
								if((solution_pre1==2)&&(solution_pre2==0))
								{
									jump_sum++;

									if(solution == 0)
									{
										badimpsnp_index[win_top]=1;
										badimpsnp_index[win_bottom]=1;
										solution = 2;
									}
									else if(solution == 2)
									{	
										if(win_top-1>=0)
											badimpsnp_index[win_top-1]=1;
										if(win_bottom+1<n_snpimp)
											badimpsnp_index[win_bottom+1]=1;
									}
								}

								if(solution==1)
								{
									if(u>=WINO)
									{
										for(int i=0;i<u;i++)
										{
											if(((hogoxn_table[0][i]=='N')||(hogoxn_table[1][i]=='N'))||((hogoxn_table[0][i]=='X')&&(hogoxn_table[1][i]=='X')))
											{
												//====================================
												//added on 03/16/2012
												imp_double_side_match++;
											//		fprintf(fpo_l,"%d\n",u);
												//====================================

												if((hogoxn_table[0][i]=='N')||(hogoxn_table[1][i]=='N'))
													imputeisxn[snpimp[indexofsnp_act[i]]]='n';
												else
													imputeisxn[snpimp[indexofsnp_act[i]]]='x';

												//====================================
												imputeisdsr[snpimp[indexofsnp_act[i]]]='d';
												//Added on 03/20/2012
												imputed_hap_0[snpimp[indexofsnp_act[i]]]=imp_hap_0;
												imputed_hap_1[snpimp[indexofsnp_act[i]]]=imp_hap_1;
												imputed_hap_total[snpimp[indexofsnp_act[i]]]=imp_hap_total;

												imputed_winsize_snp[snpimp[indexofsnp_act[i]]]=u;
												imputed_winsize_pos[snpimp[indexofsnp_act[i]]] = hgr_pos[snpimp[indexofsnp_act[u-1]]] - hgr_pos[snpimp[indexofsnp_act[0]]];
												imputed_maf_cycno[snpimp[indexofsnp_act[i]]] = impute_maf;
												imputed_maf_ref[snpimp[indexofsnp_act[i]]] = impute_maf;
												imputed_haps_sum[snpimp[indexofsnp_act[i]]] = imp_haps_sum;

												imputed_win_homo[snpimp[indexofsnp_act[i]]] = imp_win_homo;
												imputed_win_hete[snpimp[indexofsnp_act[i]]] = imp_win_hete;
												imputed_win_xx[snpimp[indexofsnp_act[i]]] = imp_win_xx;
												imputed_win_nn[snpimp[indexofsnp_act[i]]] = imp_win_nn;

												imputed_jump_sum[snpimp[indexofsnp_act[i]]] = jump_sum;
												imputed_edge_dist[snpimp[indexofsnp_act[i]]] = (u-1-i)>i?i:u-1-i;




												//====================================
												
											/*fu		printf("Imputing   %s %d   %c%c <-- %c%c\n",
													        hgr_rs[snpimp[indexofsnp_act[i]]],hgr_pos[snpimp[indexofsnp_act[i]]],hogo_xn[snpimp[indexofsnp_act[i]]][0],hogo_xn[snpimp[indexofsnp_act[i]]][1],
															imputor_refinetable[matindex_0][i],  imputor_refinetable[matindex_1][i]);
											*/	
												
												hogo_xn[snpimp[indexofsnp_act[i]]][0] = imputor_refinetable[matindex_0][i];
												hogo_xn[snpimp[indexofsnp_act[i]]][1] = imputor_refinetable[matindex_1][i];											
											}
										}
									}

									imputecycle = false;									
								}
								else if(solution == 0)
								{
									win_half = win_half-1;
									imputecycle = true;
								}
								else if(solution == 2)
								{
									win_half = win_half+2;
									imputecycle = true;
								}
								else
									printf("Error type 3");

							solution_pre1 = solution_pre2;
							solution_pre2=solution;

						}while((imputecycle==true)&&(lessimputemaxlen==true));
				}
			}
			cmphogoxn_issame = hogoxn_cmp(hogo_xn,hogo_xn_pre,snpimp,n_snpimp);
		}while(!cmphogoxn_issame);

		free(snpimp);
	}
	//==================================================================
//	fclose(fpo_l);
	//Added on 08/15/2011
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*fu		fpo=fopen("imputedhaplotype_2.txt","w");	
		fprintf(fpo,"%s", headofhaplo);
		for(int i=0;i<snpsum_hgr;i++)
		{
			if(i==snpsum_hgr-1)
				fprintf(fpo,"%s %d %c %c", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
			else
				fprintf(fpo,"%s %d %c %c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
		}
		fclose(fpo);
	*/	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//==================================================================
	//Added for rare allele impute on 01/11/2011
	int n_snpimp0=0;
	n_snpimp0=snpsum_hgr;
	int *snpimp0=(int*)malloc(snpsum_hgr*sizeof(int));
	for(int h=0;h<snpsum_hgr;h++)
		snpimp0[h]=h;

	do{
		hogoxn_cpy(hogo_xn, hogo_xn_pre,snpimp0,n_snpimp0);
		//Impute processing
		//NX locator
		for(int s=0;s<n_snpimp0;s++)
		{
			if((hogo_xn[snpimp0[s]][0]=='X')&&(hogo_xn[snpimp0[s]][1]=='X'))
			{
				jump_sum = 0;

				if(maf[snpimp0[s]]!=0)
				{
					bool imputecycle=false;
					bool lessimputemaxlen=false;
						//ini badimpsnp_index
						for(int t=0;t<n_snpimp0;t++)
							badimpsnp_index[t]=0;
					//Imputor
					//ini win=INIWIN
					int solution_pre1=0,solution_pre2=0;
					int win_half=WINI/2, win_top,win_bottom;
					do{
							win_top=s-win_half;
							win_bottom=s+win_half;
							if(win_top<0)
								win_top=0;
							if(win_bottom>n_snpimp0-1)
								win_bottom=n_snpimp0-1;

							if(win_bottom-win_top==0)
								break;
							//hogoxn_table
							int u=0;
							for(int i=0;i<2;i++)
							{
								u=0;
								for(int q=0;q<win_bottom-win_top+1;q++)
								{
									if(badimpsnp_index[win_top+q]==0)
									{
										hogoxn_table[i][u] = hogo_xn[snpimp0[win_top+q]][i];
										indexofsnp_act[u] = win_top+q;
										u++;
									}
								}
								hogoxn_table[i][u]='\0';
							}

							if(u<MAXCOMPARELEN)
								lessimputemaxlen = true;
							else
								lessimputemaxlen = false;
							
							//imputor_table for minor allele
							int r=0;
							for(int i=0;i<ref_num;i++)
							{								
								if(r_hgralleles[s][i]==main_allele[s][0]||r_hgralleles[s][i]==main_allele[s][1])
								{
									u=0;
									for(int q=0;q<win_bottom-win_top+1;q++)
									{
										if(badimpsnp_index[win_top+q]==0)
										{
											imputor_table[r][u]=r_hgralleles[snpimp0[win_top+q]][i];
											u++;
										}
									}
									imputor_table[r][u]='\0';
									r++;
								}
							}
							//====================================
							int matindex_0=0,matindex_1=0;
							int solution = match_hr0(imputor_table,r, hogoxn_table,u,imputor_refinetable,imputor_refinetable_f,hxnmatch_0,hxnmatch_1,&matindex_0,&matindex_1);	

							if((solution_pre1==2)&&(solution_pre2==0))
							{
								jump_sum++;

								if(solution == 0)
								{
									badimpsnp_index[win_top]=1;
									badimpsnp_index[win_bottom]=1;
									solution = 2;
								}
								else if(solution == 2)
								{	
									if(win_top-1>=0)
										badimpsnp_index[win_top-1]=1;
									if(win_bottom+1<n_snpimp0)
										badimpsnp_index[win_bottom+1]=1;
								}
							}

							if((solution==3)||(solution==4))
							{
							
								if(u>=WINO)
								{
									for(int i=0;i<u;i++)
										{																		
												if((hogoxn_table[0][i]=='X')&&(hogoxn_table[1][i]=='X'))
												{													
													imputeisxn[snpimp0[indexofsnp_act[i]]]='x';

													//====================================
														imputeisdsr[snpimp0[indexofsnp_act[i]]]='s';
												//Added on 03/20/2012
												imputed_hap_0[snpimp0[indexofsnp_act[i]]]=imp_hap_0;
												imputed_hap_1[snpimp0[indexofsnp_act[i]]]=imp_hap_1;
												imputed_hap_total[snpimp0[indexofsnp_act[i]]]=imp_hap_total;

												imputed_winsize_snp[snpimp0[indexofsnp_act[i]]]=u;
												imputed_winsize_pos[snpimp0[indexofsnp_act[i]]] = hgr_pos[snpimp0[indexofsnp_act[u-1]]] - hgr_pos[snpimp0[indexofsnp_act[0]]];
												imputed_maf_cycno[snpimp0[indexofsnp_act[i]]] = 0;
												imputed_maf_ref[snpimp0[indexofsnp_act[i]]] = 0;
												imputed_haps_sum[snpimp0[indexofsnp_act[i]]] = imp_haps_sum;

												imputed_win_homo[snpimp0[indexofsnp_act[i]]] = imp_win_homo;
												imputed_win_hete[snpimp0[indexofsnp_act[i]]] = imp_win_hete;
												imputed_win_xx[snpimp0[indexofsnp_act[i]]] = imp_win_xx;
												imputed_win_nn[snpimp0[indexofsnp_act[i]]] = imp_win_nn;

												imputed_jump_sum[snpimp0[indexofsnp_act[i]]] = jump_sum;
												imputed_edge_dist[snpimp0[indexofsnp_act[i]]] = (u-1-i)>i?i:u-1-i;




												//====================================

													if(solution==3)
													{
														if(imputor_refinetable[matindex_0][i]==g_hgralleles[snpimp0[indexofsnp_act[i]]][0])
														{
															//====================================
															//added on 03/16/2012
															imp_single_side_match++;
															//====================================
												/*fu			printf("Imputing   %s %d   %c%c <-- %c%c\n",
																hgr_rs[snpimp0[indexofsnp_act[i]]],hgr_pos[snpimp0[indexofsnp_act[i]]],hogo_xn[snpimp0[indexofsnp_act[i]]][0],hogo_xn[snpimp0[indexofsnp_act[i]]][1],
																g_hgralleles[snpimp0[indexofsnp_act[i]]][0], g_hgralleles[snpimp0[indexofsnp_act[i]]][1]);
*/
															hogo_xn[snpimp0[indexofsnp_act[i]]][0] =g_hgralleles[snpimp0[indexofsnp_act[i]]][0];
															hogo_xn[snpimp0[indexofsnp_act[i]]][1] = g_hgralleles[snpimp0[indexofsnp_act[i]]][1];			
														}
														else
														{
															//====================================
															//added on 03/16/2012
															imp_single_side_match++;
															//====================================
											/*fu				printf("Imputing   %s %d   %c%c <-- %c%c\n",
																hgr_rs[snpimp0[indexofsnp_act[i]]],hgr_pos[snpimp0[indexofsnp_act[i]]],hogo_xn[snpimp0[indexofsnp_act[i]]][0],hogo_xn[snpimp0[indexofsnp_act[i]]][1],
																g_hgralleles[snpimp0[indexofsnp_act[i]]][1], g_hgralleles[snpimp0[indexofsnp_act[i]]][0]);
*/
															hogo_xn[snpimp0[indexofsnp_act[i]]][0] =g_hgralleles[snpimp0[indexofsnp_act[i]]][1];
															hogo_xn[snpimp0[indexofsnp_act[i]]][1] = g_hgralleles[snpimp0[indexofsnp_act[i]]][0];		
														}
													}
													else //solution==4
													{
														if(imputor_refinetable[matindex_1][i]==g_hgralleles[snpimp0[indexofsnp_act[i]]][0])
														{
															//====================================
															//added on 03/16/2012
															imp_single_side_match++;
															//====================================
													/*fu		printf("Imputing   %s %d   %c%c <-- %c%c\n",
																hgr_rs[snpimp0[indexofsnp_act[i]]],hgr_pos[snpimp0[indexofsnp_act[i]]],hogo_xn[snpimp0[indexofsnp_act[i]]][0],hogo_xn[snpimp0[indexofsnp_act[i]]][1],
																g_hgralleles[snpimp0[indexofsnp_act[i]]][1], g_hgralleles[snpimp0[indexofsnp_act[i]]][0]);
*/
															hogo_xn[snpimp0[indexofsnp_act[i]]][0] =g_hgralleles[snpimp0[indexofsnp_act[i]]][1];
															hogo_xn[snpimp0[indexofsnp_act[i]]][1] = g_hgralleles[snpimp0[indexofsnp_act[i]]][0];			
														}
														else
														{
															//====================================
															//added on 03/16/2012
															imp_single_side_match++;
															//====================================
													/*fu		printf("Imputing   %s %d   %c%c <-- %c%c\n",
																hgr_rs[snpimp0[indexofsnp_act[i]]],hgr_pos[snpimp0[indexofsnp_act[i]]],hogo_xn[snpimp0[indexofsnp_act[i]]][0],hogo_xn[snpimp0[indexofsnp_act[i]]][1],
																g_hgralleles[snpimp0[indexofsnp_act[i]]][0], g_hgralleles[snpimp0[indexofsnp_act[i]]][1]);
*/
															hogo_xn[snpimp0[indexofsnp_act[i]]][0] =g_hgralleles[snpimp0[indexofsnp_act[i]]][0];
															hogo_xn[snpimp0[indexofsnp_act[i]]][1] = g_hgralleles[snpimp0[indexofsnp_act[i]]][1];		
														}
													}												
											}
										}
								}

								imputecycle = false;
							}
							else if(solution == 0)
							{
								win_half = win_half-1;
								imputecycle = true;
							}
							else if(solution == 2)
							{
								win_half = win_half+2;
								imputecycle = true;
							}
							else
								printf("Error type 3");

							solution_pre1 = solution_pre2;
							solution_pre2=solution;

					}while((imputecycle==true)&&(lessimputemaxlen==true));
				}
			}
		}

		cmphogoxn_issame = hogoxn_cmp(hogo_xn,hogo_xn_pre,snpimp0,snpsum_hgr);
	}while(!cmphogoxn_issame);
	free(snpimp0);

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*fu		fpo=fopen("imputedhaplotype_3.txt","w");	
		fprintf(fpo,"%s", headofhaplo);
		for(int i=0;i<snpsum_hgr;i++)
		{
			if(i==snpsum_hgr-1)
				fprintf(fpo,"%s %d %c %c", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
			else
				fprintf(fpo,"%s %d %c %c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
		}
		fclose(fpo);
fu */	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//==================================================================
		//================================================================
//fu		fpo=fopen("imputed_radom.txt","w");
//fu		fprintf(fpo,"rsNo pos\n");



		
		//================================================================
	//1> set win = 0
	for(int i=0;i<snpsum_hgr;i++)
	{
			if(((hogo_xn[i][0]=='N')||(hogo_xn[i][1]=='N'))||((hogo_xn[i][0]=='X')&&(hogo_xn[i][1]=='X')))
			{
				//Imputor
				if((hogo_xn[i][0]=='N')||(hogo_xn[i][1]=='N'))
				{
					//====================================
					//added on 03/16/2012
					imp_random_maf_cnt++;
					//====================================

					imputeisdsr[i]='r';
					imputeisxn[i]='n';
					imputed_hap_0[i] = 0;
					imputed_hap_1[i] = 0;
					imputed_hap_total[i]=ref_num;

					imputed_winsize_snp[i]=0;
					imputed_winsize_pos[i] = 0;
					imputed_maf_cycno[i] = 0;
					imputed_maf_ref[i] = 0;
					imputed_haps_sum[i] = 0;

					imputed_win_homo[i] = 0;
					imputed_win_hete[i] = 0;
					imputed_win_xx[i] = 0;
					imputed_win_nn[i] = 0;

					imputed_jump_sum[i] = 0;
					imputed_edge_dist[i] = 0;


			//fu		printf("Imputing   %s %d   %c%c <-- %c%c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1], main_allele[i][0], main_allele[i][0]);

					if(hogo_xn[i][0]=='N')
					hogo_xn[i][0] = main_allele[i][0];
					if(hogo_xn[i][1]=='N')
					hogo_xn[i][1] = main_allele[i][0];	

			//		fprintf(fpo,"%s %d 3\n", hgr_rs[i],hgr_pos[i]);
				}
				else if((hogo_xn[i][0]=='X')&&(hogo_xn[i][1]=='X'))
				{
					imputeisdsr[i]='r';
					imputeisxn[i]='x';
					imputed_hap_0[i] = 0;
					imputed_hap_1[i] = 0;
					imputed_hap_total[i]=ref_num;
					imputed_winsize_snp[i]=0;
					imputed_winsize_pos[i] = 0;
					imputed_maf_cycno[i] = 0;
					imputed_maf_ref[i] = 0;
					imputed_haps_sum[i] = 0;

					imputed_win_homo[i] = 0;
					imputed_win_hete[i] = 0;
					imputed_win_xx[i] = 0;
					imputed_win_nn[i] = 0;

					imputed_jump_sum[i] = 0;
					imputed_edge_dist[i] = 0;


			//fu		printf("Imputing   %s %d   %c%c <-- %c%c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1],g_hgralleles[i][0],g_hgralleles[i][1]) ;

					//hogo_xn[i][0] = g_hgralleles[i][1];
					//hogo_xn[i][1] = g_hgralleles[i][0];
					//or
					hogo_xn[i][0] = g_hgralleles[i][0];
					hogo_xn[i][1] = g_hgralleles[i][1];	

					if(g_hgralleles[i][0] == g_hgralleles[i][1])
					{
						//====================================
						//added on 03/16/2012
						imp_random_maf_notcnt++;
						//====================================
	//fu					fprintf(fpo,"%s %d\n", hgr_rs[i],hgr_pos[i]);
					}
					else
					{
						//====================================
						//added on 03/16/2012
						imp_random_maf_cnt++;
						//====================================
		//fu				fprintf(fpo,"%s %d\n", hgr_rs[i],hgr_pos[i]);
					}

				}
			}
	}

// fu	fclose(fpo);
	//Free g_hgralleles
	for(int i=0;i<snpsum_hgr;i++)
		free(g_hgralleles[i]);
	free(g_hgralleles);

	int d=0;
	while(headofhaplo[d] != '\n')
		d++;
	int e=0;
	//Output imputeinfo file
	
//fu	fpo=fopen("imputeinfo.txt","w");
	
	char appendstr_1[] = " another refmaf ori(1-2) isxn discrepancy\n";
	e=0;
	while(appendstr_1[e]!='\0')
	{
		headofhaplo[d+e]=appendstr_1[e];
		e++;
	}
	headofhaplo[d+e]='\0';
//fu	fprintf(fpo,"%s", headofhaplo);
/*fu	for(int i=0;i<snpsum_hgr;i++)
	{
		fprintf(fpo,"%s %d %c %c %.4f %c-%c %c %c\n", hgr_rs[i], hgr_pos[i], hogo_xn[i][0], hogo_xn[i][1], maf[i],hogo_xn_ori[i][0],hogo_xn_ori[i][1], imputeisxn[i],discrepancy[i]);
	}
	fclose(fpo);
*/
	//Output imputedhaplotype file and preimputehaplotype file
	char appendstr_2[]=" another\n";
	e=0;
	while(appendstr_2[e]!='\0')
	{
		headofhaplo[d+e]=appendstr_2[e];
		e++;
	}
	headofhaplo[d+e]='\0';

//	fpo=fopen("imputedhaplotype.txt","w");	
	fpo=fopen(impute,"w");	
	fprintf(fpo,"%s", headofhaplo);
	for(int i=0;i<snpsum_hgr;i++)
	{
		if(i==snpsum_hgr-1)
			fprintf(fpo,"%s %d %c %c", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
		else
			fprintf(fpo,"%s %d %c %c\n", hgr_rs[i],hgr_pos[i],hogo_xn[i][0],hogo_xn[i][1]);
	}
	fclose(fpo);

/*fu	fpo=fopen("preimputehaplotype.txt","w");	
	fprintf(fpo,"%s", headofhaplo);
	for(int i=0;i<snpsum_hgr;i++)
	{
		if(i==snpsum_hgr-1)
			fprintf(fpo,"%s %d %c %c", hgr_rs[i],hgr_pos[i],hogo_xn_ori[i][0],hogo_xn_ori[i][1]);
		else
			fprintf(fpo,"%s %d %c %c\n", hgr_rs[i],hgr_pos[i],hogo_xn_ori[i][0],hogo_xn_ori[i][1]);
	}
	fclose(fpo);
	//==========================================================================
	//Added on 03/20/2012
	fpo=fopen("infosum.txt","w");	
	fprintf(fpo,"rsNo pos xn dsr fre_h0 fre_h1 fre_total winsize_snp winsize_pos maf_cycno maf_ref haps_sum win_homo win_hete jump_sum edge_dist win_xx win_nn\n");
	for(int i=0;i<snpsum_hgr;i++)
	{
		fprintf(fpo,"%s %d %c %c %d %d %d %d %d %.2f %.2f %d %d %d %d %d %d %d\n",hgr_rs[i],hgr_pos[i],imputeisxn[i],imputeisdsr[i],imputed_hap_0[i],imputed_hap_1[i],imputed_hap_total[i],
						imputed_winsize_snp[i],imputed_winsize_pos[i],imputed_maf_cycno[i],imputed_maf_ref[i],imputed_haps_sum[i],imputed_win_homo[i],imputed_win_hete[i],imputed_jump_sum[i],imputed_edge_dist[i],imputed_win_xx[i],imputed_win_nn[i]);
		\

	}
	fclose(fpo);
	//==========================================================================
	//Added on 03/23/2012
 	fpo=fopen("quality_score.txt","w");	
 */
	fpo=fopen(qscore,"w");
	fprintf(fpo,"rsNo pos QS_geno QS_haplo\n");
	for(int i=0;i<snpsum_hgr;i++)
	{
		getQScore(imputeisxn[i],imputeisdsr[i],imputed_hap_0[i],imputed_hap_1[i],imputed_hap_total[i],imputed_winsize_snp[i],imputed_winsize_pos[i],imputed_maf_cycno[i],imputed_haps_sum[i],imputed_win_homo[i],imputed_win_hete[i],
											imputed_jump_sum[i],imputed_edge_dist[i],imputed_win_xx[i]);
		fprintf(fpo,"%s %d %.2f %.2f\n",hgr_rs[i],hgr_pos[i],q_score_geno,q_score_haplo);
	}
	fclose(fpo);



	//==========================================================================
	//Added on 03/16/2012
	int imp_total = imp_double_side_match+imp_single_side_match+imp_random_maf_cnt+imp_random_maf_notcnt;
/* fu
	fpo=fopen("match_stype.txt","w");
	fprintf(fpo,"No solution: %d\n",imp_random_maf_cnt+imp_random_maf_notcnt);
	fprintf(fpo,"No single solution: %d\n",imp_single_side_match+(imp_random_maf_cnt+imp_random_maf_notcnt));
	fprintf(fpo,"No big gap: %d\n",imp_total-imp_double_side_match-imp_single_side_match-(imp_random_maf_cnt+imp_random_maf_notcnt));
	fprintf(fpo,"Total imputed: %d\n",imp_total);
	fprintf(fpo,"Ref_homo_imputeXX: %d\n",imp_random_maf_notcnt);
	fprintf(fpo,"Ref_nothomo_imputeXX_or_imputeNN: %d\n",imp_random_maf_cnt);
	fclose(fpo);
*/
	//Free r_hgralleles
	for(int i=0;i<snpsum_hgr;i++)
		free(r_hgralleles[i]);
	free(r_hgralleles);
	//Free imputor_table
	for(int i=0;i<ref_num;i++)
		free(imputor_table[i]);
	free(imputor_table);
	//Free imputor_refinetable
	for(int i=0;i<ref_num;i++)
		free(imputor_refinetable[i]);
	free(imputor_refinetable);
	//Free imputor_refinetable_f
	for(int i=0;i<ref_num;i++)
		free(imputor_refinetable_f[i]);
	free(imputor_refinetable_f);
	//Free hogoxn_table
	for(int i=0;i<2;i++)
		free(hogoxn_table[i]);
	free(hogoxn_table);
	//Allocate hogo_xn_ori
	for(int i=0;i<snpsum_hgr;i++)
		free(hogo_xn_ori[i]);
	free(hogo_xn_ori);
	//Allocate hogo_xn
	for(int i=0;i<snpsum_hgr;i++)
		free(hogo_xn[i]);
	free(hogo_xn);
	//Free badimpsnp_index
	free(badimpsnp_index);
	//Free imputeisxn
	free(imputeisxn);
	//Free discrepancy
	free(discrepancy);
	//Free indexofsnp_act
	free(indexofsnp_act);
	//Free hxnmatch_0
	free(hxnmatch_0);
	free(hxnmatch_1);

	//Free hgr_pos
	free(hgr_pos);
	//Free hgr_rs
	for(int i=0;i<snpsum_hgr;i++)
		free(hgr_rs[i]);
	free(hgr_rs);
	//Free h_index
	free(h_index);
	//Free main_allele
	for(int i=0;i<snpsum_hgr;i++)
		free(main_allele[i]);
	free(main_allele);
	//Free maf
	free(maf);
	//Free hogo_xn_pre
	for(int i=0;i<snpsum_hgr;i++)
		free(hogo_xn_pre[i]);
	free(hogo_xn_pre);

	//==========================================================================
	//Added on 03/20/2012
	//1)
	free(imputed_hap_0);
	free(imputed_hap_1);
	free(imputed_hap_total);
	//2)
	free(imputed_winsize_snp);
	//3)
	free(imputed_winsize_pos);
	//4)
	free(imputed_maf_cycno);
	//5)
	free(imputed_maf_ref);
	//6)
	free(imputed_haps_sum);
	//7)
	free(imputed_win_homo);
	free(imputed_win_hete);	
	//8)
	free(imputed_jump_sum);
	//9)
	free(imputed_edge_dist);
	//10)
	free(imputed_win_xx);
	free(imputed_win_nn);

	//==========================================================================

	printf("*************done!*************\n");	
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;
	printf("Runing time: %.2f seconds\n", duration);
	printf("*******************************\n");
	//Output running file
/*fu	FILE *fpr;
	fpr=fopen("runingtime.txt","w");	
	fprintf(fpr,"Runing time: %.2f seconds\n",duration);
	fclose(fpr);
*/
	system("pause");
	return 0;
}

int *array_union(int *a, int na, int *b, int nb, int *nc)
{
	  int *c=NULL;
	  int i=0,j=0,t=0;;

       while (i<na&&j<nb)
       {
		    if (a[i]<b[j])
			{
				 if (!c) 
					 c=(int*)malloc(sizeof(int));
				 else 
					 c=(int*)realloc(c,sizeof(int)*(t+1));
				 c[t++]=a[i++];
			}
			 else if (b[j]<a[i]) 
			 {
				if (!c) 
					c=(int*)malloc(sizeof(int));
				 else 
					 c=(int*)realloc(c,sizeof(int)*(t+1));
				 c[t++]=b[j++];
			}
			else if (a[i]==b[j]) 
			{
				 if (i<j) i++;
				 else j++;
			  }
	   }

	    while (i<na)
		  {
			  c=(int*)realloc(c,sizeof(int)*(t+1));
			  c[t++]=a[i++];
		  }
		  while (j<nb)
		   {
				 c=(int*)realloc(c,sizeof(int)*(t+1));
				 c[t++]=b[j++];
		   }
		   *nc=t;
		   return c;	
}

int cmp_rise(const void *a, const void *b)
{
	return *(char*)a-*(char*)b;
}

bool hogoxn_cmp(char **a, char **b, int *c, int nc)
{
	for(int i=0;i<nc;i++)
	{
		if((a[c[i]][0]!=b[c[i]][0])||(a[c[i]][1]!=b[c[i]][1]))
			return false;
	}
	return true;
}

void hogoxn_cpy(char **a, char **b, int *c, int nc)
{
	for(int i=0;i<nc;i++)
	{
		b[c[i]][0] = a[c[i]][0];
		b[c[i]][1] = a[c[i]][1];		
	}
}
int match_hr0(char **hr, int n_hr, char **hxn, int l_hxn, char **hr_refine, char** hr_refine_f, int *hxnmat_0, int *hxnmat_1,int *matind_0, int *matind_1)
{
	imp_hap_total = n_hr;
	getHxnInf(hxn,l_hxn);

	qsort(hr, n_hr, sizeof(char*),cmporder);

	strcpy(hr_refine[0],hr[0]);
	int j=1;
	for(int i=0;i<n_hr-1;i++)
	{
		if(strcmp(hr[i], hr[i+1])!=0)
		{
			strcpy(hr_refine[j],hr[i+1]);
			j++;
		}
	}

	imp_haps_sum = j;

	//For hxnmatch_0	
	cpytable_hr(hr_refine_f,hr_refine,j);
	setxntohr(hr_refine_f,hxn[0],j,l_hxn);
	int hmatsum_0 = cmphxnhr(hr_refine_f,hxn[0],j,hxnmat_0);
	//For hxnmatch_1	
	cpytable_hr(hr_refine_f,hr_refine,j);
	setxntohr(hr_refine_f,hxn[1],j,l_hxn);
	int hmatsum_1 = cmphxnhr(hr_refine_f,hxn[1],j,hxnmat_1);

	//compare match
	int matchnum = 0;
	int r_t = 0, v_t=0;

	if((hmatsum_0!=0)&&(hmatsum_1!=0))
		matchnum=2;
	else if((hmatsum_0==0)&&(hmatsum_1==0))
		matchnum=0;
	else if((hmatsum_0!=0)&&(hmatsum_1==0))
	{
		//=========================================================
		//Added on 03/20/2012
	//	imp_hap_0 = getHapFre(hr_refine[r_t],hr,n_hr);
	//	imp_hap_1 = getHapFre(hr_refine[v_t],hr,n_hr);		
		//=========================================================
		if(hmatsum_0==1)
		{
			matchnum=3;
			r_t=hxnmat_0[0];

			imp_hap_0 = getHapFre(hr_refine[r_t],hr,n_hr);
			imp_hap_1 = 0;
		}
		else
			matchnum=2;
	}
	else
	{
		if(hmatsum_1==1)
		{
			matchnum=4;
			v_t=hxnmat_1[0];

			imp_hap_1 = getHapFre(hr_refine[v_t],hr,n_hr);		
			imp_hap_0 = 0;
		}
		else
			matchnum=2;
	}

	*matind_0=r_t;
	*matind_1=v_t;

	return matchnum;
}
int match_hr(char **hr, int n_hr, char **hxn, int l_hxn, char **hr_refine, char** hr_refine_f, int *hxnmat_0, int *hxnmat_1,int *matind_0, int *matind_1)
{
	imp_hap_total = n_hr;
	getHxnInf(hxn,l_hxn);

	qsort(hr, n_hr, sizeof(char*),cmporder);

	strcpy(hr_refine[0],hr[0]);
	int j=1;
	for(int i=0;i<n_hr-1;i++)
	{
		if(strcmp(hr[i], hr[i+1])!=0)
		{
			strcpy(hr_refine[j],hr[i+1]);
			j++;
		}
	}

	imp_haps_sum = j;

	//For hxnmatch_0	
	cpytable_hr(hr_refine_f,hr_refine,j);
	setxntohr(hr_refine_f,hxn[0],j,l_hxn);
	int hmatsum_0 = cmphxnhr(hr_refine_f,hxn[0],j,hxnmat_0);
	//For hxnmatch_1	
	cpytable_hr(hr_refine_f,hr_refine,j);
	setxntohr(hr_refine_f,hxn[1],j,l_hxn);
	int hmatsum_1 = cmphxnhr(hr_refine_f,hxn[1],j,hxnmat_1);

	//compare match
	int matchnum = 0;
	int r_t = 0, v_t=0;
	for(int r=0;r<hmatsum_0;r++)
		for(int v=0;v<hmatsum_1;v++)
		{
			bool match = matchtest(hr_refine[hxnmat_0[r]],hr_refine[hxnmat_1[v]],hxn, l_hxn);
			if(match == true)
			{
				matchnum++;
				r_t = hxnmat_0[r];
				v_t= hxnmat_1[v];
			}		
		}
		*matind_0=r_t;
		*matind_1=v_t;
		
		if(matchnum >2)
			matchnum = 2;
		//=========================================================
		//Added on 03/20/2012
		if(matchnum == 1)
		{
			imp_hap_0 = getHapFre(hr_refine[r_t],hr,n_hr);
			imp_hap_1 = getHapFre(hr_refine[v_t],hr,n_hr);
		}


		//=========================================================
	return matchnum;
}
bool matchtest(char *a_1, char *a_2, char **b, int l)
{
	bool matchresult = true;
	for(int i=0;i<l;i++)
	{
		if((b[0][i]=='X')&&(b[1][i]=='X'))
		{
			if(a_1[i]==a_2[i])
			{
				matchresult = false;
				break;
			}
		}
	}
 return matchresult;
}
int cmphxnhr(char **a, char *b, int j, int *c)
{
	int matchnum=0;
	for(int i=0;i<j;i++)
	{
		if(strcmp(b, a[i])==0)
		{
			c[matchnum]=i;
			matchnum++;
		}
	}
	return matchnum;
}
void setxntohr(char **a, char *b,int j, int k)
{
	for(int i=0;i<k;i++)
	{
		if((b[i]=='X')||(b[i]=='N'))
		{
			for(int m=0;m<j;m++)
				a[m][i]=b[i];
		}
	}
}
int cmporder(const void *a, const void *b)
{
	return strcmp(*(char**)a, *(char**)b);
}
void cpytable_hr(char **a, char **b, int j)
{
	for(int i=0;i<j;i++)
	{
		strcpy(a[i],b[i]);
	}
}

int getHapFre(char *a, char**b,int len)
{
	int sum = 0;

	for(int i=0;i<len;i++)
	{
		if(strcmp(a,b[i]) == 0)
			sum++;
	}

	return sum;
}

void getHxnInf(char** a, int len)
{
	int homo=0;
	int hete = 0;
	int xx = 0;
	int nn = 0;

	for(int i=0;i<len;i++)
	{
		if((a[0][i]=='X'||a[0][i]=='x')&&(a[1][i]=='X'||a[1][i]=='x'))
			xx++;
		else if((a[0][i]=='N'||a[0][i]=='n')&&(a[1][i]=='N'||a[1][i]=='n'))
			nn++;
		else if(a[0][i]==a[1][i])
			homo++;
		else if(a[0][i]!=a[1][i])
			hete++;
	}

	imp_win_homo = homo;
	imp_win_hete = hete;
	imp_win_xx = xx;
	imp_win_nn = nn;
}

void getQScore(char xn, char dsr, int fre_h0, int fre_h1, int fre_total, int winsize_snp, int winsize_pos, double maf_cycno, int haps_sum,int win_homo, int win_hete,int jump_sum, int edge_dist, int win_xx)
{
	double s_haplo = 0;
	double	s_geno = 0;
	double t_subscore_freq_h = 0;
	double t_winsize_snp = 0;
	double t_winsize_pos = 0;
	double t_maf_cycno = 0;
	double t_haps_sum = 0;
	double t_win_homo = 0;
	double t_win_hete = 0;
	double t_jump_sum = 0;
	double t_edge_dist = 0;
	double t_win_xx = 0;

	//======================================
	//1> t_subscore_freq_h
	int freq_sum = fre_h0+fre_h1;
	int freq_min = fre_h0 < fre_h1? fre_h0 : fre_h1;
	if(freq_sum == 0)
		t_subscore_freq_h = 0;
	else
	{
		if(freq_min == 0)
		{
			if(freq_sum >= 13)
				t_subscore_freq_h = 1.5;
			else
				t_subscore_freq_h = (double)(freq_sum-1)/12.0 + 0.5;
		}
		else if(freq_min == 1)
		{
			if(freq_sum==2)
				t_subscore_freq_h=1.9;
			else if(freq_sum>=3 && freq_sum<=12)
				t_subscore_freq_h=(double)(freq_sum-3)/10.0 + 2.0;
			else if(freq_sum>= 13)
				t_subscore_freq_h = 3.0;
		}
		else if(freq_min >= 2)
		{
			if(freq_sum>=4 && freq_sum<=14)
				t_subscore_freq_h = ((double)(freq_sum-4)/11.0)*1.9 + 3.1;
			else if(freq_sum >= 15)
				t_subscore_freq_h = 5.0;
		}
	}

	//======================================
	//2> t_winsize_snp
	if(winsize_snp==10 || winsize_snp==11)
		t_winsize_snp = -5.0;
	else if(winsize_snp==12 || winsize_snp==13)
		t_winsize_snp = -1.0;
	else if(winsize_snp==14 || winsize_snp==15)
			t_winsize_snp = 3.0;
	else if(winsize_snp==16 || winsize_snp==17)
			t_winsize_snp = 5.0;
	else if(winsize_snp==18 || winsize_snp==19)
			t_winsize_snp = 3.0;
	else if(winsize_snp >= 20)
		   t_winsize_snp = 4.0;

	//======================================
	//3> t_winsize_pos
	if(winsize_pos >= 300000)
		t_winsize_pos = 0;
	else if(winsize_pos>=100001 && winsize_pos<=299999)
		t_winsize_pos = 4.0 - (((double)(winsize_pos-100001)/200000)*3.0);
	else if(winsize_pos>=0 && winsize_pos<=100000)
		t_winsize_pos = 5.0;
	//======================================
	//4> t_maf_cycno
	if(maf_cycno>=0 && maf_cycno<=0.14)
		t_maf_cycno = (maf_cycno/0.14)*0.40;
	else if(maf_cycno>=0.15 && maf_cycno<=0.24)
		t_maf_cycno = 5.0;
	else if(maf_cycno>=0.25 && maf_cycno<=0.50)
		t_maf_cycno = 0.4-((maf_cycno-0.25)/0.26)*0.24;
	//======================================
	//5> t_haps_sum
	int haps_each = (int)( (double)fre_total/(double)haps_sum);
	if(haps_each <= 2)
		t_haps_sum = 0;
	else if(haps_each >= 20)
		t_haps_sum = 5;
	else
		t_haps_sum = ((double)(haps_each-2)/18.0)*3.0 + 2.0;

	//======================================
	//6> t_win_homo
	if(win_homo < (winsize_snp/5))
		t_win_homo = 5;
	else
		t_win_homo = 5.0*(double)(winsize_snp/5)/(double)win_homo;

	//======================================
	//7> t_win_hete
	if(win_hete == 0)
		t_win_hete = 5.0;
	else if(win_hete==1)
		t_win_hete = -10.0;
	else if(win_hete==2)
		t_win_hete = 1.0;
	else if(win_hete>=3 && win_hete<=7)
		t_win_hete = ((double)(win_hete-3)/5.0)*3.0 + 2.0;
	else if(win_hete >= 8)
		t_win_hete = 5.0;

		//======================================
	//8> t_jump_sum
	if(jump_sum ==0)
		t_jump_sum = 5.0;
	else if(jump_sum ==1)
		t_jump_sum = 4.0;
	else if(jump_sum ==2)
		t_jump_sum = 3.0;
	else if(jump_sum >=3 && jump_sum <=7)
		t_jump_sum = 2.5-((double)(jump_sum-3)/5.0)*2;
	else if(jump_sum >= 8)
		t_jump_sum = 0;
	//======================================
	//8> t_edge_dist
	if(edge_dist < winsize_snp/4)
		t_edge_dist = 5.0;
	else
		t_edge_dist = 5.0*(double)(winsize_snp/4)/(double)edge_dist;
	//======================================
	//8> t_win_xx
	if(win_xx > winsize_snp/4)
		t_win_xx = 5;
	else
		t_win_xx = 5.0*(double)win_xx/(double)(winsize_snp/4);


	if(xn == 'x')
	{
		s_haplo = (t_subscore_freq_h*0.37 + t_winsize_snp*0.1 + t_haps_sum*0.37 + t_win_homo*0.03 + t_win_hete*0.1 + t_jump_sum*0.02 + t_win_xx*0.01)/(5.0);
		s_geno = 1.00;
	}
	else if(xn=='n')
	{
		s_haplo = (t_subscore_freq_h*0.37 + t_winsize_snp*0.1 + t_haps_sum*0.37 + t_win_homo*0.03 + t_win_hete*0.1 + t_jump_sum*0.02 + t_win_xx*0.01)/(5.0);
		s_geno = (t_subscore_freq_h*0.24 + t_winsize_snp*0.12 + t_winsize_pos*0.065 + t_maf_cycno*0.065 + t_haps_sum*0.25 +  t_win_homo*0.004 +  t_win_hete*0.07 + t_jump_sum*0.19 + t_edge_dist*0.003 + t_win_xx*0.003)/(5.0);
	}
	else
	{
		s_haplo = 1.00;
		s_geno = 1.00;
	}

	if(dsr == 'r')
	{
		if(xn == 'x')
		{
			s_haplo = 0;
			s_geno = 1.0;
		}
		else if(xn == 'n')
		{
			s_haplo = 0;
			s_geno = 0.8;
		}
	}
	
	q_score_haplo = s_haplo;
	q_score_geno = s_geno;
}
