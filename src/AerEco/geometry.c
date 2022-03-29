#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (void)
{
   int i,j,k,g,i1,j1,t,smod = 0,s = 0,sum = 0;
   int kobj,kmod,sign,nfm,i11,i22;
   double x,r,e;
//------------------------------------Calculation of the number of frames-------------------------------------
printf("begin programm\n");
   char *test = "geodat.dat";
   FILE *fin = fopen(test,"r");
   fscanf(fin,"%d",&kobj);
   for (i = 0; i < kobj; i++) {
      fscanf(fin,"%d",&kmod);
      smod = smod + kmod;
      for (j = 0; j < kmod; j++) {
         fscanf(fin,"%d",&sign);
         fscanf(fin,"%d",&nfm);
         sum = sum + nfm;
         fscanf(fin,"%d %d",&i11,&i22);
         nfm = nfm*12;
         for (k = 0; k < nfm; k++) {
            fscanf(fin,"%lf",&x);
         }
      }
   }
   fclose(fin);
//--------------------------------------------Reading constants-----------------------------------------------
printf("begin read options\n");
   char *opt = "Options.txt";
   fin = fopen(opt,"r");
   fscanf(fin,"%lf",&e);
   fclose(fin);
//------------------------------------------List of various points--------------------------------------------
printf("begin read points\n");
   int be[2][smod];
   double root[3][4][sum];
   sum = 0;
   t = 0;
   fin = fopen(test,"r");
   fscanf(fin,"%d",&kobj);
   for (i = 0; i < kobj; i++) {
      fscanf(fin,"%d",&kmod);
      for (j = 0; j < kmod; j++) {
         fscanf(fin,"%d",&sign);
         fscanf(fin,"%d",&nfm);
         fscanf(fin,"%d %d",&i11,&i22);
         be[0][t] = sum;
         for (k = 0; k < nfm; k++) {
            for (g = 0; g < 4; g++) {
               fscanf(fin,"%lf %lf %lf",&root[0][g][sum],&root[1][g][sum],&root[2][g][sum]);                  //Simple reading
            }
            sum++;          //number of frames
         }
         be[1][t] = sum - 1;
         t++;
      }
   }
   fclose(fin);
   
   double dr1,dr2,dr0;
   dr1 = 10000000.0;
   dr2 = 10000000.0;
   for(i = 0; i < sum; i++) {
       for (j = 0; j < 4; j++) {
	  g = j+1;
	  if(g == 4) {g = 0;}
	  dr0 = ( (root[0][j][i] - root[0][g][i])*(root[0][j][i] - root[0][g][i]) + (root[1][j][i] - root[1][g][i])*(root[1][j][i] - root[1][g][i]) +  (root[2][j][i] - root[2][g][i])*(root[2][j][i] - root[2][g][i]));
	  dr0 = sqrt(dr0);
	  if(dr0 != 0.0) {
	       if(dr0 < dr1){dr1 = dr0;}
	       if(dr0 < dr2 && dr0 > e) {dr2 = dr0;}
	  }
       }	 
   }     
printf(" distanse between adjasent points, unequal points  %lf %lf \n",dr1,dr2);
   
   
   
printf("begin formig unuq point list \n");
char *pnt = "Points.txt";
FILE *fout = fopen(pnt,"w");
   s = 0;
   double point[3][sum*4];
   for (i = 0; i < sum; i++) {
      for (j = 0; j < 4; j++) {
         t = 0;
         for (i1 = 0; i1 <= i; i1++) {
            for (j1 = 0; j1 < 4; j1++) {
               r = (root[0][j][i] - root[0][j1][i1])*(root[0][j][i] - root[0][j1][i1]) + (root[1][j][i] - root[1][j1][i1])*(root[1][j][i] - root[1][j1][i1]) + (root[2][j][i] - root[2][j1][i1])*(root[2][j][i] - root[2][j1][i1]);
               if (r < e*e) {
                  if ( (i1 < i) || ((i1 == i) && (j1 < j)) ) {
                     t = 1;
                     break;
                  }
               }
            }
            if (t == 1) {
               break;
            }
         }
         if (t == 0) {
fprintf(fout,"%lf %lf %lf\n",root[0][j][i],root[1][j][i],root[2][j][i]);
            point[0][s] = root[0][j][i];
            point[1][s] = root[1][j][i];
            point[2][s] = root[2][j][i];
            s++;              //number of different points
         }
      }
   }
fclose(fout);
printf("Number of different points is: s = %d,\nTotal number of frames is: sum = %d\n",s,sum);

//----------------------------------------------Numbered list of cells----------------------------------------
   int frame[4][sum];
char *NumCell = "NumCell.txt";
fout = fopen(NumCell,"w");
   for (i = 0; i < sum; i++) {
      for (j = 0; j < 4; j++) {
         for (k = 0; k < s; k++) {
           r = ( (root[0][j][i] - point[0][k])*(root[0][j][i] - point[0][k]) + (root[1][j][i] - point[1][k])*(root[1][j][i] - point[1][k]) + (root[2][j][i] - point[2][k])*(root[2][j][i] - point[2][k]) );
	   if ( r < e*e ) {
fprintf(fout,"%d ",k + 1);
               frame[j][i] = k + 1;
               break;
            }
         }
      }
fprintf(fout,"\n");
   }
fclose(fout);
//----------------------------------------------List of different ribs----------------------------------------
   int rib[2][sum*4];
   int t0,t1,t2,t3,nseg;
   j = 0;
   for (i = 0; i < sum; i++) {
      t0 = 0;
      t1 = 0;
      t2 = 0;
      t3 = 0;
      for (k = 0; k < j; k++) {
         if ( ((frame[0][i] == rib[0][k]) && (frame[1][i] == rib[1][k])) || ((frame[1][i] == rib[0][k]) && (frame[0][i] == rib[1][k])) ) {
            t0 = 1;
         }
         if ( ((frame[1][i] == rib[0][k]) && (frame[2][i] == rib[1][k])) || ((frame[2][i] == rib[0][k]) && (frame[1][i] == rib[1][k])) ) {
            t1 = 1;
         }
         if ( ((frame[2][i] == rib[0][k]) && (frame[3][i] == rib[1][k])) || ((frame[3][i] == rib[0][k]) && (frame[2][i] == rib[1][k])) ) {
            t2 = 1;
         }
         if ( ((frame[3][i] == rib[0][k]) && (frame[0][i] == rib[1][k])) || ((frame[0][i] == rib[0][k]) && (frame[3][i] == rib[1][k])) ) {
            t3 = 1;
         }
      }

      if (t0 == 0) {
	if (frame[0][i] != frame[1][i]) {
           rib[0][j] = frame[0][i];
           rib[1][j] = frame[1][i];
           j++;
	}
      }
      if (t1 == 0) {
	if (frame[1][i] != frame[2][i]) {
           rib[0][j] = frame[1][i];
           rib[1][j] = frame[2][i];
           j++;
	}
      }
      if (t2 == 0) {
	if (frame[2][i] != frame[3][i]) {
           rib[0][j] = frame[2][i];
           rib[1][j] = frame[3][i];
           j++;
	}
      }
      if (t3 == 0) {
	if (frame[3][i] != frame[0][i]) {
           rib[0][j] = frame[3][i];
           rib[1][j] = frame[0][i];
           j++;
	}
      }
   }
   nseg = j;
printf("Successfully, number of different segments is: nseg = %d\n",nseg);
//----------------------------------------------In what frame do the ribs lie---------------------------------
   int ribcell[5][nseg];
   for (i = 0; i < nseg; i++) {
      ribcell[0][i] = -1;
      ribcell[1][i] = -1;
      ribcell[2][i] = -1;
      ribcell[3][i] = -1;
      ribcell[4][i] = -1;
  }

   for (i = 0; i < sum; i++) {
      frame[0][i] = frame[0][i];
      frame[1][i] = frame[1][i];
      frame[2][i] = frame[2][i];
      frame[3][i] = frame[3][i];
      for (j = 0; j < nseg; j++) {
         if (  ((rib[0][j] == frame[0][i]) && (rib[1][j] == frame[1][i])) || ((rib[0][j] == frame[1][i]) && (rib[1][j] == frame[2][i])) || ((rib[0][j] == frame[2][i]) && (rib[1][j] == frame[3][i])) || ((rib[0][j] == frame[3][i]) && (rib[1][j] == frame[0][i])) || ((rib[0][j] == frame[1][i]) && (rib[1][j] == frame[0][i])) || ((rib[0][j] == frame[2][i]) && (rib[1][j] == frame[1][i])) || ((rib[0][j] == frame[3][i]) && (rib[1][j] == frame[2][i])) || ((rib[0][j] == frame[0][i]) && (rib[1][j] == frame[3][i])) ) {
            ribcell[4][j] = ribcell[3][j];
            ribcell[3][j] = ribcell[2][j];
            ribcell[2][j] = ribcell[1][j];
            ribcell[1][j] = ribcell[0][j];
            ribcell[0][j] = i;
         }
      }
   }
   char *adj = "Adjacent_0.txt";
   fout = fopen(adj,"w");
   for(i=0;i<nseg;i++)   {
         fprintf(fout,"%d %d %d %d %d \n",ribcell[0][i]+1,ribcell[1][i]+1,ribcell[2][i]+1,ribcell[3][i]+1,ribcell[4][i]+1);
   }
   fclose(fout);
   
printf("Successfully filled array ribcell[2][nseg]\n");
//--------------------------------------Corresponding numbering-----------------------------------------------
   char *af = "AdjacentFrames.txt";
   int r1,r2,r3,r4,r5,p1,p2,js1,is1;
   fout = fopen(af,"w");
   is1 = 0;
   js1 = 0;

   for (i = 0; i < nseg; i++) 
   {
      r1 = ribcell[0][i];
      r2 = ribcell[1][i];
      r3 = ribcell[2][i];
      r4 = ribcell[3][i];
      r5 = ribcell[4][i];
      fprintf(fout,"%d %d %d %d %d ",r1 + 1,r2 + 1,r3+1,r4+1,r5+1);

       
      p1 = rib[0][i];
      p2 = rib[1][i];
      if(r2 == -1) {js1++;}
      if(r3 == -1 && r2 != -1 ) {is1++;}
      if(r4 == -1 && r3 != -1 ) {is1 = is1 + 2;}
      if(r5 == -1 && r4 != -1 ) {is1 = is1 + 3;}
      if(r5 != -1 ) {is1 = is1 + 4;}
 
      if (r1 != -1) 
      {
         for (j = 0; j < 4; j++) 
	 {
            if (p1 == frame[j][r1]) 
	    {
               fprintf(fout,"%d ",j + 1);
;
	       break;
            }
         }
         for (j = 0; j < 4; j++) 
	 {
            if (p2 == frame[j][r1]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
      } 
      else 
      {
         fprintf(fout," 0 0");
      }

      if (r2 != -1) 
      {
         for (j = 0; j < 4; j++) 
	 {
            if (p1 == frame[j][r2]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
         for (j = 0; j < 4; j++) 
	 {
            if (p2 == frame[j][r2]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
      } 
      else 
      {
         fprintf(fout," 0 0");
      }
 

       if (r3 != -1) 
      {
         for (j = 0; j < 4; j++) 
	 {
            if (p1 == frame[j][r3]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
         for (j = 0; j < 4; j++) 
	 {
            if (p2 == frame[j][r3]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
      } 
      else 
      {
         fprintf(fout," 0 0");
      }
      
      if (r4 != -1) 
      {
         for (j = 0; j < 4; j++) 
	 {
            if (p1 == frame[j][r4]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
         for (j = 0; j < 4; j++) 
	 {
            if (p2 == frame[j][r4]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
      } 
      else 
      {
         fprintf(fout," 0 0");
      }
      
      if (r5 != -1) 
      {
         for (j = 0; j < 4; j++) 
	 {
            if (p1 == frame[j][r5]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
         for (j = 0; j < 4; j++) 
	 {
            if (p2 == frame[j][r5]) 
	    {
               fprintf(fout,"%d ",j + 1);
               break;
            }
         }
      } 
      else 
      {
         fprintf(fout," 0 0");
      }

      fprintf(fout,"\n");
   }
   fclose(fout);
   
// doubled segments print   
   int n1,n2,n3,n4,n5,k1,k2,k3,k4,k5;
   af = "AdjacentFrames.txt";
   FILE * f1 = fopen(af,"r");   
   char *af1 = "Adjacent_Doubled_Segments.txt";
   FILE * f2 = fopen(af1,"w");   
     for (i = 0; i < nseg; i++) {
          fscanf(f1,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&r1,&r2,&r3,&r4,&r5,&n1,&k1,&n2,&k2,&n3,&k3,&n4,&k4,&n5,&k5);	
	  if(r3 != 0) {
	     fprintf(f2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",r1,r2,r3,r4,r5,n1,k1,n2,k2,n3,k3,n4,k4,n5,k5);	
	  }	    
     }       
   
   
  fclose(f1);
  fclose(f2);
   
//--------------------------------------------------RESULT----------------------------------------------------
   char *res = "data.rwg";
   FILE *fres = fopen(res,"w");


   
   fprintf(fres,"%d\n",s);
   for (i = 0; i < s; i++) {
      fprintf(fres,"%lf %lf %lf\n",point[0][i],point[1][i],point[2][i]);
   }

   fprintf(fres,"%d\n",sum);
   for (i = 0; i < sum; i++) {
      fprintf(fres,"%d %d %d %d\n",frame[0][i],frame[1][i],frame[2][i],frame[3][i]);
   }

   fprintf(fres,"%d\n",is1);
   fin = fopen(af,"r");
   for (i = 0; i < nseg; i++) {
      fscanf(fin,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&r1,&r2,&r3,&r4,&r5,&n1,&k1,&n2,&k2,&n3,&k3,&n4,&k4,&n5,&k5);
      if ( r2 != 0 ) {
         fprintf(fres,"%d %d %d %d %d %d\n",r1,r2,n1,k1,n2,k2);
      }
      if ( r3 != 0 ) {
         fprintf(fres,"%d %d %d %d %d %d\n",r1,r3,n1,k1,n3,k3);
      }
      if ( r4 != 0 ) {
         fprintf(fres,"%d %d %d %d %d %d\n",r1,r4,n1,k1,n4,k4);
      }
      if ( r5 != 0 ) {
         fprintf(fres,"%d %d %d %d %d %d\n",r1,r5,n1,k1,n5,k5);
      }

   }
   fclose(fin);

   fprintf(fres,"%d\n",js1);
   fin = fopen(af,"r");
   for (i = 0; i < nseg; i++) {
      fscanf(fin,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&r1,&r2,&r3,&r4,&r5,&n1,&k1,&n2,&k2,&n3,&k3,&n4,&k4,&n5,&k5);
      if ((r2 == 0) ) {
         fprintf(fres,"%d %d %d %d %d %d\n",r1,r2,n1,k1,n2,k2);
      }
   }
   fclose(fin);

   fprintf(fres,"%d\n",smod);
   for (i = 0; i < smod; i++) {
      fprintf(fres,"%d %d\n",be[0][i] + 1,be[1][i] + 1);
   }

   fclose(fres);
   return 0;
}
