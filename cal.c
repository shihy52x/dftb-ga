#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>

#define TARGET          27
#define T               300
#define pi              3.1415926
#define Max_Atom         200
#define CMax_Atom        500

typedef struct
{
 double  fitness;
 double  gen[Max_Atom][3];
 double  cgen[CMax_Atom][3];
 double  ep;
}ga_struct;


void rand_init(int natm,double *x, double boxl)

{
double RR=0.0;
int i=0;
int ii=0;
int j=0;
double Rfcc=1.7;
double rmin=0.7*Rfcc;
boxl=pow(natm,0.333)*2.7*1.5;
  for ( i=0;i<natm;i++)
    {
      for ( j=0;j<3;j++)
        {
          *(x+i*3+j)=(rand()/(RAND_MAX*1.0)*2-1)*boxl/2.0;
        }

      for( ii=0;ii<i;ii++)
        {
          RR=0.0;
          RR=   pow(*(x+i*3+0)-*(x+(ii)*3+0),2);
          RR=RR+pow(*(x+i*3+1)-*(x+(ii)*3+1),2);
          RR=RR+pow(*(x+i*3+2)-*(x+(ii)*3+2),2);
          RR=sqrt(RR);
//          printf ("i=%d,ii=%d,RR=%f\n",i,ii,RR);
          if (RR<rmin )
            {i=i-1;
//          printf ("i=%d,ii=%d,RR=%f\n",i,ii,RR);
            }
         }
   }
}




int cmp(const void *a,const void *b)
{
      return ((double *)b)[0] >((double*)a)[0] ? 1:-1;
}


//Function: center: centerize the geometryand rotate around the center

void
center(int natm ,double a[Max_Atom][3])
{
 double x0=0,y0=0,z0=0;
 double xx=0,yy=0,zz=0,xxold=0;
// random rotation angle generator, rotate about y axis and z axis
//
 double  dtheta1=rand()%360;
 double  dtheta2=rand()%360;
         dtheta1=dtheta1/360.0*2*pi;
         dtheta2=dtheta2/360.0*2*pi;
         



 int i=0;
// calculation of center coordination x0,y0,z0
  for(i=0;i<natm;i++)
    {
      x0=a[i][0]+x0;
      y0=a[i][1]+y0;
      z0=a[i][2]+z0;
    }
  x0=x0/natm;
  y0=y0/natm;
  z0=z0/natm;

// translation of center to 0,0,0 and rotation
//
  for(i=0;i<natm;i++)
    {   
      xx=a[i][0]-x0;
      yy=a[i][1]-y0;
      zz=a[i][2]-z0;
      xxold=xx;
      xx=xx*cos(dtheta1)-yy*sin(dtheta1);
      yy=xxold*sin(dtheta1)+yy*cos(dtheta1);
      xxold=xx;
      xx=xx*cos(dtheta2)-zz*sin(dtheta2);
      zz=xxold*sin(dtheta2)+zz*cos(dtheta2);
      a[i][0]=xx;
      a[i][1]=yy;
      a[i][2]=zz;
//    printf("after  rotation 1  %12.4f %12.4f %12.4f\n",a[i][0],a[i][1],a[i][2]);
    } 

//Function: mate:

}


void
shift(int natm ,double a[Max_Atom][3],double ddptc)
{
 double zmin=0;
 int i=0;
//initialize geometry center 

  for(i=0;i<natm;i++)
   {
   
   if(zmin>a[i][2]) {zmin=a[i][2];} 
//   printf("before center 1  %12.4f %12.4f %12.4f\n",a[i][0],a[i][1],a[i][2]);
   }
 
double  dgap=ddptc-(zmin-0); //0 is the z coord of graphne sheet 
   
// translation of cluster 
//
  for(i=0;i<natm;i++)
   {
    a[i][2]=a[i][2]+dgap;//dgap is the atuall distance. ddptc is required distance
   }

}


void
init_carbon_i(int cnatm,int POPSIZE, ga_struct *beta_population,int i)
{

double *cx = (double *) malloc(3*cnatm*sizeof(double));
int j=0;
FILE  *fp_input=fopen("c.coord","r");
read_coord(fp_input,cnatm,cx);
  for (j=0;j<cnatm;j++)
    {
      beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
      beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
      beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
    }
free (cx);
fclose(fp_input);
  
}


void
init_carbon(int cnatm,int POPSIZE, ga_struct *beta_population)
{
double *cx = (double *) malloc(3*cnatm*sizeof(double));
int i=0,j=0;
FILE  *fp_input=fopen("c.coord","r");
read_coord(fp_input,cnatm,cx);
for(i=0;i<POPSIZE;i++)
  {
    for (j=0;j<cnatm;j++)
      {
        beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
        beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
        beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
  
      }
  }
free(cx);
fclose(fp_input);
  
}



void
init_population(int ptnatm, int cnatm, int POPSIZE,ga_struct *population, ga_struct *beta_population,double boxl,int flag_res)
{
 double *x = (double *) malloc(3*ptnatm*sizeof(double));
 double *cx = (double *) malloc(3*cnatm*sizeof(double));
 double efinal=0;
 FILE  *fp_coordc=fopen("c.coord","r");
 FILE  *fp_coordpt=fopen("pt.coord","r");

 if (flag_res==1&&fp_coordpt==NULL)
    {printf("ERROR, can not open pt.coord!!\n");exit(0);}
//read coordinates of graphene from c.coord

 read_coord(fp_coordc,cnatm,cx);
 int i=0;
 int j=0;


 for(i=0;i<POPSIZE;i++)
   {
     population[i].fitness=0;
     beta_population[i].fitness=0;

//reading from first generation
//
//read_coord(fp_input2,natm,x);
     if(flag_res==0)
        {rand_init(ptnatm,x,boxl);}
     else if(flag_res==1)
        {
          read_coord(fp_coordpt,ptnatm,x);
        }


     for (j=0;j<ptnatm;j++)
       {
//                printf("%d %lf %lf %lf \n", skip,x,y,z);
		
         population[i].gen[j][0]=(double)(x[j*3+0]);
         population[i].gen[j][1]=(double)(x[j*3+1]);
         population[i].gen[j][2]=(double)(x[j*3+2]);
         beta_population[i].gen[j][0]=0;
         beta_population[i].gen[j][1]=0;
         beta_population[i].gen[j][2]=0;

        }

      for (j=0;j<cnatm;j++)
        {
          population[i].cgen[j][0]=(double)(cx[j*3+0]);
          population[i].cgen[j][1]=(double)(cx[j*3+1]);
          population[i].cgen[j][2]=(double)(cx[j*3+2]);
          beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
          beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
          beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
         }

//       center(natm,population[i].gen);  
          printf("%d %lf %lf %lf \n",i,population[i].cgen[cnatm-1][0],population[i].cgen[cnatm-1][1],population[i].cgen[cnatm-1][2]);
          printf("%d %lf %lf %lf \n",i,population[i].gen[0][0],population[i].gen[ptnatm-1][1],population[i].gen[ptnatm-1][2]);
     }
    free (x);
    free (cx);
    fclose(fp_coordc);
    if(fp_coordpt!=NULL){ fclose(fp_coordpt);}

}



double grep_energy(){
   char line[1024];
   double energy;
   FILE *fp =fopen("detailed.out","r");
   int n=1;
   while(1){
    if (fgets(line,1024,fp) ==NULL) n=0;
    if(n==0) break;
    if(strstr(line,"SCC is NOT converged")) printf("Not Converged! "); 
    if(strstr(line,"Total energy")){
      char *pch;
      pch=strtok(line," ");
      int i=0;
      while (pch !=NULL){
         i=i+1;
         if(i==5) energy=strtod(pch,NULL);
       pch=strtok (NULL," ");
      }
    }
  } 

   fclose(fp);
   return energy;                             
}

double cal_energy(double *xx, double *yy, int n1, int n2, int nc, double *cell){
    double energy;
    int status;
    FILE *fp2 ,*fp3;
    fp2=fopen("geo.gen","w"); write_gen(fp2,n1,n2,nc,xx,yy,cell,"Pt","Ru","C"); fclose(fp2);
    status=system("dftb+>dftb_screen");
    energy=grep_energy();
    fp3 =fopen("geo_end.gen","r"); read_gen(fp3,n1,n2,nc,xx,yy);fclose(fp3); 
    printf("Energy= %12.4f\n",energy);
    return grep_energy();                             
}




void
cal_pop_energy(int POPSIZE,ga_struct *population,int ptnatom, int runatom, int cnatom, double *cell ,int argc, char **argv)
{
  double efinal=0;
  int i=0;

  for(i=0;i<POPSIZE;i++)
    { 
//       cal_energy(&efinal,population[i].gen,population[i].cgen,0,0,ptnatom, cnatom ,argc, argv);       
       efinal=cal_energy(population[i].gen,population[i].cgen,ptnatom,runatom,cnatom,cell);
       population[i].ep=efinal;
    }
 
}





void
normal_fitness(int POPSIZE,ga_struct *population) 
{      

double  mine=population[0].ep;
double  maxe=population[POPSIZE-1].ep;
int     i=0;

    for (i=0;i<POPSIZE;i++)
      {        
        if(mine==maxe)
          {
         population[i].fitness=1.0;
          } 
        else
          {
          population[i].fitness=1-0.7*(population[i].ep-mine)/(maxe-mine);
        
          }
      }

}
 

int
sort_func (const void *e1,const void *e2)

{
	return (((ga_struct *) e1)->ep > ((ga_struct *)e2)->ep ? 1:-1);
}



// function count: count the number of atoms below or above xx plane
int
count_x(int natm, double a[Max_Atom][3],int xx)
{

int  i=0;

  for(i=0;i<natm;i++)
   {
	    if( a[i][0] >0)
		 {xx=xx+1;   
 		 }
   }
 return xx;
}

//Function: mate:


double min(double a,double b)
{
if(a>b){return b;}else{return a;}
}

void 
elitism(int esize,int ptnatm, int cnatm ,ga_struct *population,ga_struct *beta_population)
{
 int i=0,j=0;
 for (i=0;i<esize;i++)
   {
//      center(natm,population[i].gen);  
      beta_population[i].ep=population[i].ep;
      for(j=0;j<ptnatm;j++)
        {
          beta_population[i].gen[j][0]=population[i].gen[j][0];
          beta_population[i].gen[j][1]=population[i].gen[j][1];
          beta_population[i].gen[j][2]=population[i].gen[j][2];
        }
      for(j=0;j<cnatm;j++)
        {
          beta_population[i].cgen[j][0]=population[i].cgen[j][0];
          beta_population[i].cgen[j][1]=population[i].cgen[j][1];
          beta_population[i].cgen[j][2]=population[i].cgen[j][2];
         }
   }
}

void 
mutate ( int natm, double Mu,int POPSIZE,ga_struct *population)
{

 int rand_pop=0;
 int rand_atom=0;
 int rand_direct=0;
 int rand_N=0;
 double rand_l=0.0;
 int i=0,j=0;
 for(i=0;i<POPSIZE*Mu;i++)
   {

     rand_pop=rand()%POPSIZE;
     rand_atom=rand()%natm;
     rand_N=rand()%50;
     for(j=0;j<rand_N;j++)
       {
         rand_l=(rand()/(RAND_MAX*1.0)-0.5)*4;
         rand_direct=rand()%3;
  //    printf("randpop= %d, rand_atom= %d rand_N= %d, rand_direct= %d rand_l= %lf\n", rand_pop,rand_atom,rand_N,rand_direct,rand_l);
          population[rand_pop].gen[rand_atom][rand_direct]=population[rand_pop].gen[rand_atom][rand_direct]+rand_l;

        }  
    }

}

void 
mutate_perm (int ptnatm, int runatm, double Mu,int POPSIZE,ga_struct *population){
 int    Nmetal=ptnatm+runatm;
 double randnum;
 int randN1,randN2;
 double temp[3];
 int i ,j,ii;
 for(i=1;i<POPSIZE;i++){ // do not mutate the 0th cluster 
   randnum=rand()/(RAND_MAX*1.0); //generage random number in [0,1]
   randN1=rand()%(ptnatm); //generage random number in [0,natm)-type1
   randN2=rand()%(runatm)+ptnatm; //generage random number in [ptnatm,Nmetal)-type2
   printf(" i=%d %lf %d %d\n",i,randnum,randN1,randN2);
   if(randnum<Mu){
     for (ii=0;ii<13;ii++) printf(" %d %lf\n",ii,population[i].gen[ii][0]);
     memcpy(temp,population[i].gen[randN1],sizeof(temp));
     memcpy(population[i].gen[randN1],population[i].gen[randN2],sizeof(temp));
     memcpy(population[i].gen[randN2],temp,sizeof(temp));
     for (ii=0;ii<13;ii++) printf(" %d %lf\n",ii,population[i].gen[ii][0]);
   }
 }

}



void 
mate( int ptnatm, int runatm, int cnatm, double *cell, int esize,int POPSIZE,ga_struct *population,ga_struct *beta_population, int Temp, double e_mate, double ddptc,int min_step,int argc, char **argv)
{
 int Nmetal=ptnatm+runatm;
 int    rand_index=0;
 int    randa=0;
 int    randb=0;
 int    min_index=0;
 double rand_p=0.0;

 int    Na=0,Nb=0;
 double plandelta=0;
 double bulkd=1.0; //the distance between two bulk after join;
 double parentrate=0.2;
 double Ne_mate=e_mate*ptnatm;
 double planex1,planex2;
 double p2copy[Nmetal][3],p1copy[Nmetal][3];
 double p2copy_2[Nmetal][3],p1copy_2[Nmetal][3];
 planex1=0;
 planex2=0;
 int    i,j,k,ii,jj;
 int    type1[Nmetal],type2[Nmetal];
 int    type1_s[Nmetal],type2_s[Nmetal];
 int    type_new[Nmetal];
 int Npt,Ndiff,type_change;
 int rank;
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 for (ii=0;ii<ptnatm ;ii++) type1[ii]=1;
 for (ii=ptnatm;ii<Nmetal;ii++) type1[ii]=2;
   FILE  *fbad=fopen("bad.xyz","a+");
// FILE  *fbad=fopen("bad.xyz","w");
// all graphene in this generation will be reset to c.coord.
 for(i=0;i<1;i++)
   {
   
// select two parents 
     for (j=0;j<2;j++)
       {	
         do
           {
             rand_index=rand()%POPSIZE;
             rand_p=rand()/(RAND_MAX*1.0);
           } while(population[rand_index].fitness<rand_p);
         if(j==0){randa=rand_index;}
         if(j==1){randb=rand_index;}
              //   if(randa==randb){j=j-1;}
        }  

     memcpy( p1copy,population[randa].gen,sizeof(p1copy));
     memcpy( p2copy,population[randb].gen,sizeof(p2copy));
     memcpy( type2,type1,sizeof(type1));
     center(Nmetal,p1copy);
     center(Nmetal,p2copy);
     memcpy( p1copy_2,p1copy,sizeof(p1copy));
     memcpy( p2copy_2,p2copy,sizeof(p2copy));
///     write_coord(fbad,ptnatm,runatm,cnatm,p1copy,beta_population[0].cgen,0);
//    write_coord(fbad,ptnatm,runatm,cnatm,p1copy,beta_population[0].cgen,0);
//     write_coord(fbad,ptnatm,runatm,cnatm,p2copy,beta_population[0].cgen,0);
     qsort(p1copy,Nmetal,sizeof(p1copy[0]),cmp);
     qsort(p2copy,Nmetal,sizeof(p2copy[0]),cmp);
///     write_coord(fbad,ptnatm,runatm,cnatm,p1copy,beta_population[0].cgen,0);
///     write_coord(fbad,ptnatm,runatm,cnatm,p2copy,beta_population[0].cgen,0);
     for (ii=0;ii<=Nmetal;ii++){
        for (j=0;j<=Nmetal;j++){
        if (p1copy[ii][0]==p1copy_2[j][0]) type1_s[ii]=type1[j];
        if (p2copy[ii][0]==p2copy_2[j][0]) type2_s[ii]=type2[j];
        }
     } 
//     for (ii=0;ii<13;ii++) printf("type1  atom[%d]= %d %lf\n",ii,type1_s[ii],p1copy[ii][0]);
     Na=count_x(Nmetal,p1copy,0);
     Nb=Nmetal-Na;
     planex1=p1copy[Na-1][0];
     planex2=p2copy[Na][0];
     plandelta=planex1-planex2;
     Npt=0;
     Ndiff=0;
     for (j=0;j<Nmetal;j++)
       {
         if(j<Na)
           { 
                beta_population[i].gen[j][0]=p1copy[j][0];
                beta_population[i].gen[j][1]=p1copy[j][1];
                beta_population[i].gen[j][2]=p1copy[j][2];
                type_new[j]=type1_s[j];
                if(type_new[j]==1) Npt=Npt+1;
           } 
         else
           {
                beta_population[i].gen[j][0]=p2copy[j][0]+plandelta-bulkd*plandelta/fabs(plandelta);
                beta_population[i].gen[j][1]=p2copy[j][1];
                beta_population[i].gen[j][2]=p2copy[j][2];
                type_new[j]=type2_s[j];
                if(type_new[j]==1) Npt=Npt+1;
            }
        }
     Ndiff=ptnatm-Npt;
     if(Ndiff>0)  type_change=2;
     if(Ndiff<0) {type_change=1; Ndiff=-Ndiff;}
     for (j=Na;j<Nmetal;j++){
         if (Ndiff==0) break;
         if (type_new[j]==type_change){
            type_new[j]=type_change%2+1;
            Ndiff--;              
          }
     } 
     
     k=0;
     ii=0;
     jj=ptnatm;
     for (j=0;j<Nmetal;j++){
//             printf("j=%d %d\n",j,type_new[j]);
             if(type_new[j]==1) {k=ii; ii++;}
             if(type_new[j]==2) {k=jj; jj++;}
                p1copy[k][0]=beta_population[i].gen[j][0];
                p1copy[k][1]=beta_population[i].gen[j][1];
                p1copy[k][2]=beta_population[i].gen[j][2];
////                printf("k=%d %lf\n",k,p1copy[j][0]);
     }
//     for (ii=0;ii<13;ii++) printf("type  atom[%d]= %d %lf\n",ii,type_new[ii],p1copy[ii][0]);

     init_carbon_i(cnatm,POPSIZE, beta_population,i);
     shift(ptnatm,beta_population[i].gen,ddptc);
     memcpy( beta_population[i].gen,p1copy,sizeof(p1copy));
     beta_population[i].ep=cal_energy(beta_population[i].gen,beta_population[i].cgen,ptnatm,runatm,cnatm,cell);
     write_coord(fbad,ptnatm,runatm,cnatm,beta_population[i].gen,beta_population[0].cgen,0);
     fflush(fbad);
 
     if(population[randa].ep<population[randb].ep)
       {min_index=randa;}
     else
       {min_index=randb;}
          	    
  
     if(beta_population[i].ep-Ne_mate<population[min_index].ep)
       {
         printf("%dth core Eold= %12.4lf  Enew= %12.4lf %d %d SUCCESSFULLY MATED !!! \n", rank,population[min_index].ep,beta_population[i].ep,randa,randb);
       }
      else  
        {
          printf("%dth core Eold= %12.4lf  Enew =%12.4f %d %d NOT          MATED !!!\n" ,rank,population[min_index].ep,beta_population[i].ep,randa,randb);
          i=i-1;          
   
        }
 }
/// fclose(fbad);
}


 






void
swap ( ga_struct ** p1, ga_struct* * p2)
{
  ga_struct *tmp = *p1;
  *p1 = *p2;
  *p2 = tmp;
}


 void filter(int natm,int POPSIZE, int *  NEWPOPSIZE, double delte,ga_struct *population)
{
  int i,j,k;
 *NEWPOPSIZE=POPSIZE;
  for (i=1;i<POPSIZE;i++)
      {
       for (j=i+1;j<POPSIZE;j++)
        {
         if(fabs(population[j].ep-population[i].ep)<delte)
           
           {
             population[j].ep=population[0].ep;
               for (k=0;k<natm;k++)
                {
                population[j].gen[k][0]=population[0].gen[k][0];
                population[j].gen[k][1]=population[0].gen[k][1];
                population[j].gen[k][2]=population[0].gen[k][2];
                }
                  
//              population[j].ep=1000;
//           *NEWPOPSIZE=*NEWPOPSIZE-1;
           
           }
       } 
      }

   qsort (population, POPSIZE, sizeof(ga_struct),sort_func);

}












//****************************************************************************************//
//             Main Function Block                                                        // 
//****************************************************************************************//


int  main(int argc, char **argv){

    int rc,num_processors,rank;
    rc=MPI_Init(&argc,&argv);
    if (rc != MPI_SUCCESS) {
    	printf("Error starting MPI program. Terminating\n");
           MPI_Abort(MPI_COMM_WORLD, rc);
    }
   
    MPI_Comm_size(MPI_COMM_WORLD,&num_processors);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   //**************************************************************************************//
   //                  reading input information by all cores                              //
   //**************************************************************************************//
    FILE     *fpenergy=fopen("./energy","w"); // Print the energy evolution
    FILE     *fprestart=fopen("./restart","w"); //print the generation and popsize coordinates
    FILE     *fpoptim=fopen("./optim.xyz","wb");
    FILE     *fp=fopen("data.txt","r");
    time_t    current_time;
    double    seconds,start,end,seconds_new,seconds_old,seconds_total;
    struct    timespec now,tmstart;
    char*     c_time_string,c_time_final;
    int       flag_converge=0;
    int       i=0; 
   
    srand(time(NULL)+rank);
    clock_gettime(CLOCK_REALTIME, &tmstart);
    printf("hello from processor %ld\n",rank);
    start = (double)clock() /(double) CLOCKS_PER_SEC;
    seconds_old=  (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9);
    current_time = time(NULL);
    c_time_string = ctime(&current_time);
    printf("Current time is %s\n", c_time_string);
    printf("hello from processor %ld\n",rank);
   
 
    char skip[10];
    FILE     *fp_input2 = fopen("ga_dftb_input1.1","r");
    int       NSTEP=0;
    int       step=0;
    double    ee_mate=0 ;       // the energy cretiria for accepting children cluster. the larger, the less restrict. can be 0.1 0 or -0.1 
    double    ELITRATE=0.2;
    int       POPSIZE=0;
    int       min_step=0;
    int       NEWPOPSIZE=POPSIZE;
    double    delte=0.00001 ; //0.00001;
    int       ptnatm,runatm,cnatm;
    double    boxl=0;
    double    Mu=0.2;
    double    dptc;
    int       glob=0,globconvg=20;
    double    globe[30];
    int       Temp; 
    int       flag_res=0;
    double    cell[9];
    int       esize=0;
    int       num_cores_child;
      
    fscanf(fp_input2,"%d  %d %d %s\n", &ptnatm, &runatm,&cnatm,&skip);
    fscanf(fp_input2,"%d %s\n ", &POPSIZE,skip);
    fscanf(fp_input2,"%d %s\n", &NSTEP,skip);
    fscanf(fp_input2,"%d %s\n", &globconvg,skip);
    fscanf(fp_input2,"%lf %s\n",&ELITRATE,skip);
    fscanf(fp_input2,"%lf %s\n",&delte,skip);
    fscanf(fp_input2,"%d %s\n",&Temp,skip);
    fscanf(fp_input2,"%d %s\n",&min_step,skip);
    fscanf(fp_input2,"%lf %s\n",&ee_mate,skip);
    fscanf(fp_input2,"%lf %s\n",&dptc,skip);
    fscanf(fp_input2,"%lf %s\n",&boxl,skip);
    fscanf(fp_input2,"%d %s\n",&flag_res,skip);
    fscanf(fp_input2,"%d %s\n",&num_cores_child,skip);
    int ii=0;
    for (ii=0;ii<3;ii++){
        fscanf(fp_input2,"%lf %lf %lf\n", cell+ii*3+0,cell+ii*3+1,cell+ii*3+2);


   //**************************************************************************************//
   //         Above are the information that all processors need to know                   //
   //**************************************************************************************//

    if(rank==0){
        }
       
        printf("********************JOB started*****************************\n");
        printf("********************JOB started*****************************\n");
        printf("********************JOB started*****************************\n\n\n");


        printf("Attention!! The  dftb excutable file must be in ~/bin and must be named 'dftb+' !!\n");
        printf("Attention!! The  dftb excutable file must be in ~/bin and must be named 'dftb+' !!\n");
        printf("Attention!! The  dftb excutable file must be in ~/bin and must be named 'dftb+' !!\n");
        printf("Attention!! The  dftb excutable file must be in ~/bin and must be named 'dftb+' !!\n");
        printf("Attention!! The  dftb excutable file must be in ~/bin and must be named 'dftb+' !!\n");


        printf("Number of atoms %d %d %d\n", ptnatm,runatm,cnatm);
        if(ptnatm+runatm>=Max_Atom || cnatm>=CMax_Atom) exitp("Number of atoms exceed allowed Max!");
        printf("POPSIZE         %d\n",POPSIZE);
        printf("NSTEP           %d\n",NSTEP);
        printf("globconvg       %d\n",globconvg);
        printf("ELITRATE        %lf\n",ELITRATE);
        printf("delte           %lf\n",delte);
        printf("Temprature      %d\n",Temp);
        printf("minimiz   step  %d\n",min_step);
        printf("ee_mate         %lf\n",ee_mate);
        printf("dptc            %lf\n",dptc);
        printf("intial boxl     %lf\n",boxl);
        printf("reading from pt_coord?  %d\n",flag_res);
        printf("Total number of processsors required  %d\n",num_processors);
        printf("Number of processors for each child  %d\n",num_cores_child);
        printf("\n\n\n******** end reading input information ***********************\n\n\n");
    }
    if(POPSIZE!=num_processors/num_cores_child){
        printf("Error!,POPSZIE=%d num_processrs=%d num_cores_child=%d num_processrs/num_cores_child=%d",\
                POPSIZE,num_processors,num_cores_child,num_processors/num_cores_child);
        MPI_Abort(MPI_COMM_WORLD,0);
    }    

 
    ga_struct *population = malloc(sizeof(ga_struct)*POPSIZE);
    ga_struct *beta_population = malloc(sizeof(ga_struct)*POPSIZE);

 
//****************************************************************************************//
//define new mpi structure for data transfer between processors                           //
//****************************************************************************************//
    int count=4;
    int length[4]={1,3*Max_Atom,3*CMax_Atom,1};
    MPI_Aint offset[4]={offsetof(ga_struct,fitness),offsetof(ga_struct,gen),offsetof(ga_struct,cgen),offsetof(ga_struct,ep)} ;
    MPI_Datatype type[4]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Datatype popDatatype;
    MPI_Type_struct(count,length,offset,type,&popDatatype);
    MPI_Type_commit(&popDatatype);
//****************************************************************************************//


//****************************************************************************************//
//Initialization the population for all processors                                        // 
//****************************************************************************************//

    
    
    if(rank==0){
        printf("Total number of processors is %d\n",num_processors);
        printf("Total number of population is %d\n",POPSIZE);
        printf ("For each candidate, there are %d processors to be used\n",num_processors/POPSIZE);
        init_population(ptnatm+runatm,cnatm,POPSIZE,population,beta_population,boxl,flag_res);
        printf("argc=%d\n",argc);
//        int  i=0,j=0;
//        int  esize=POPSIZE*ELITRATE;
        cal_pop_energy(POPSIZE,population,ptnatm,runatm,cnatm,cell,argc, argv);
///        for (i=0;i<POPSIZE;i++){
////////      center(ptnatm+runatm,population[i].gen);  
///           write_coord(fpoptim,ptnatm,runatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
//           shift(ptnatm+runatm, population[i].gen,dptc);  
///        } 
//        for (i=0;i<POPSIZE;i++)
//           write_coord(fpoptim,ptnatm,runatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
        fflush(fpoptim);
    }

    char command[30];
    char filename[30];
 
    sprintf(command,"mkdir core%d",rank);
    system(command);
    sprintf(command,"cp dftb_in.hsd *.coord *.skf core%d",rank);
    system(command);
    sprintf(filename,"core%d",rank);
    chdir(filename);
    system("pwd");    

//****************************************************************************************//
//              Loop started                                                              // 
//****************************************************************************************//

    for (step=0;step<NSTEP;step++){ 


       //*********************************************************************************//
       //              master core preparation                                            // 
       //*********************************************************************************//

        if(rank==0){
            cal_pop_energy(POPSIZE,population,ptnatm,runatm,cnatm,cell, argc, argv);
            printf("\n\n\n\n***********************************************\n");
            printf(  "***********************************************\n");
            printf(  "***********************************************\n");
            printf("Gen= %d starting optimization..................\n\n",step); 
        
            qsort (population, POPSIZE, sizeof(ga_struct),sort_func);
        
            normal_fitness(POPSIZE,population);
        
            for(i=0;i<POPSIZE;i++){
                fprintf(fpenergy,"%d num %d   %lf\n",step, i,population[i].ep);
                fflush(fpenergy);
            }  
        
            printf("fabs %lf\n",fabs(population[0].ep-population[POPSIZE-1].ep));
        
            if(fabs(population[0].ep-population[POPSIZE-1].ep)<delte){    
                fprintf(fpenergy,"%d %lf  %lf global minimum \n",step,population[0].ep,population[0].fitness);
                globe[glob]=population[POPSIZE-1].ep;
                glob=glob+1;           
            }
        
            if(glob>20){
                if(fabs(globe[glob]-globe[glob-20])<delte){
                    fprintf(fpenergy,"%d %lf  %lf final global minimum \n",step,population[0].ep,population[0].fitness);
                    flag_converge=1;
                }
            }

///// PPreserve the first esize parentes from the previous generation. copy from population to beta_population
            esize=POPSIZE*ELITRATE;
            if (esize<1){esize=1;}
            elitism(esize,ptnatm+runatm, cnatm,population,beta_population);

            for (i=1;i<num_processors;i++)
                MPI_Send(population,POPSIZE,popDatatype,i,0,MPI_COMM_WORLD); 
        
            
               //send coordinates and energy information to other ith  processors
        }else MPI_Recv(population,POPSIZE,popDatatype,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

       //*********************************************************************************//
       //              master core preparation ended                                      // 
       //*********************************************************************************//



       //*********************************************************************************//
       //              other cores mating start                                           // 
       //*********************************************************************************//


       //receive coordinates and energy information from 0(master) core
        MPI_Bcast(&flag_converge, 1, MPI_INT,0,MPI_COMM_WORLD);
//        printf("rank=%d,xyz=%lf\n",rank,population[0].gen[0][0]);

        if(flag_converge ==1) {MPI_Finalize(); exit(0);}    
        
//        for (i=0;i<POPSIZE;i++)
//           write_coord(fpoptim,ptnatm,runatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
        fflush(fpoptim);
//        printf("rank=%d,xyz=%lf\n",rank,population[0].gen[0][0]);

       // Generate the rest part of beta_generation by mating process
        if(rank!=0&&rank<POPSIZE)
            mate(ptnatm,runatm,cnatm,cell,esize,POPSIZE,population,beta_population,Temp,ee_mate,dptc,min_step,argc, argv);
        MPI_Barrier(MPI_COMM_WORLD);
                       
        if(rank!=0&&rank<POPSIZE)
            MPI_Send(&beta_population[0],1,popDatatype,0,1,MPI_COMM_WORLD); //send the information of first children(optimized) from other processors to master core 

       //*********************************************************************************//
       //              other cores mating ended                                           // 
       //*********************************************************************************//
 


        if(rank==0){
            for (i=1;i<POPSIZE;i++)
                MPI_Recv(&population[i],1,popDatatype,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //recieve coordinates and energy information to other ith  processors
            if (ptnatm!=0&&runatm!=0)
                mutate_perm (ptnatm,runatm, Mu, POPSIZE, beta_population);
            fprestart=fopen("./restart","w");
            for (i=0;i<POPSIZE;i++){
                write_coord(fpoptim,ptnatm,runatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
                write_coord(fprestart,ptnatm,runatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
                fflush(fpoptim);
                fflush(fprestart);
            } 
            fclose(fprestart);

         
        clock_gettime(CLOCK_REALTIME, &now);
        seconds_new = (double)((now.tv_sec+now.tv_nsec*1e-9));
        seconds=seconds_new-seconds_old;
        seconds_total = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
        printf("\nWall time for this generation is %lf s\n",seconds);
        printf("\nWall time totally  %lf s\n",seconds_total);
        seconds_old=seconds_new;
    
        }
 	  
   }




//*********************************************************************************
//************************************loop ended***********************************
//*********************************************************************************


    if(rank==0){ 
        printf("\n********************JOB FINISHED*****************************\n");
        fclose(fpenergy);
        fclose(fpoptim);
       
///////// time information
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        printf("Current time is %s\n", c_time_string);
        
         // measure elapsed wall time
        clock_gettime(CLOCK_REALTIME, &now);
        seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
        printf("wall time %fs\n", seconds);
       
         // measure CPU  time
        end = (double)clock() / (double) CLOCKS_PER_SEC;
        printf("cpu time %fs\n", end - start);
        printf("\n********************JOB FINISHED*****************************\n");
        free(population);
        free(beta_population);
    }     
}
