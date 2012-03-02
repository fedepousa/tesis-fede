#define Federico using
#define Javier namespace
#define Pousa std;
#include <ilcplex/cplex.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#define PI 3.14159265358979323846





Federico Javier Pousa


static int CPXPUBLIC mycutcallback(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);

static void free_and_null (char** ptr);

struct Punto{
	double x;
	double y;
	Punto(){}
	Punto(double a, double b){
		x = a;
		y = b;
	}
};

double norm(Punto &p1, Punto &p2){
	return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
Punto resta(Punto &p1, Punto &p2){
	return Punto(p1.x-p2.x,p1.y-p2.y);
}

double prod(Punto &p1, Punto &p2){
	return p1.x*p2.x+p1.y*p2.y;
}

int main(int argc, char **argv){

	//~ alfa es el peso correspondiente a tsp, beta a angtsp
	double alfa = 0.5;
	double beta = 0.5;
	int solstat;
	int status;
	double objval;
	double *x = NULL;
	double *y = NULL;
	double *pi = NULL;
	double *slack = NULL;
	double *dj = NULL;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int i, j;
	int cur_numrows, cur_numcols;
	int cantPuntos;
	
	
	//~ if(argc != 2){
		//~ uso();
		//~ goto TERMINATE;
	//~ }



///////////////////////////////////////////////////////////////////////////////////////////////
//Lectura de la instancia
	FILE *fin;
	fin = fopen("entrada", "r");
	fscanf(fin, "%d", &cantPuntos);
	//~ Muestra la cantidad de puntos
	//~ cerr << cantPuntos << endl;
	#define MAXPUNTOS 100
	Punto puntos[MAXPUNTOS];
	for(int i=0;i<cantPuntos;i++){
		fscanf(fin, "%lf %lf", &puntos[i].x, &puntos[i].y);
	}
	
	//~ Muestra los puntos levantados
	//~ for(int i=0;i<cantPuntos;i++){
		//~ cerr << i << " " << puntos[i].x << " " << puntos[i].y << endl;
	//~ }
	int indice = 0;
	//Este array seria para saber de que variable estoy hablando
	string revind[(MAXPUNTOS*MAXPUNTOS*MAXPUNTOS)<<1];
	int indx[MAXPUNTOS][MAXPUNTOS];
	double dist[MAXPUNTOS][MAXPUNTOS];
	for(int i=0;i<cantPuntos;i++){
		for(int j=0;j<cantPuntos;j++){
			dist[i][j] = norm(puntos[i], puntos[j]);
			stringstream ss;
			ss << i << " " << j;
			revind[indice]=ss.str();
			indx[i][j] = indice++;
			//Muestra las distancias
			//~ cerr << i << " " << j << " " << dist[i][j] << endl;
		}
	}
	
	int indy[MAXPUNTOS][MAXPUNTOS][MAXPUNTOS];
	double ang[MAXPUNTOS][MAXPUNTOS][MAXPUNTOS];
	for(int i=0;i<cantPuntos;i++){
		for(int j=0;j<cantPuntos;j++){
			for(int k=0;k<cantPuntos;k++){
				ang[i][j][k] = 0;
				stringstream ss;
				ss << i << " " << j << " " << k;
				revind[indice]=ss.str();
				indy[i][j][k] = indice++;
			}
		}
	}
	Punto aux1, aux2;
	for(int i=0;i<cantPuntos;i++){
		for(int j=0;j<cantPuntos;j++){
			if(i==j)continue;
			for(int k=0;k<cantPuntos;k++){
				//~ if(j==k||i==k)continue;
				if(k==j)continue;
				if(k==i){
					ang[i][j][k]=PI;
					continue;
				}
				aux1 = resta(puntos[j], puntos[i]);
				aux2 = resta(puntos[k], puntos[j]);
				
				ang[i][j][k] = acos(prod(aux1,aux2)/(norm(puntos[i], puntos[j])*norm(puntos[j],puntos[k])));
				//~ Muestra en grados:
				//~ cerr << i << " " << j << " " << k << " " << ang[i][j][k]*180.0/PI << endl;
			}
		}
	}
	fclose(fin);
	
	
	
////////////////////////////////////////////////////////////////////////////////////////////////
//Para mostrar como quedan las matrices de indices
	//~ for(int i=0;i<cantPuntos;i++){
		//~ for(int j=0;j<cantPuntos;j++){
			//~ cerr << indx[i][j] << " ";
		//~ }
		//~ cerr << endl;
	//~ }
	
	//~ for(int i=0;i<cantPuntos;i++){
		//~ for(int j=0;j<cantPuntos;j++){
			//~ for(int k=0;k<cantPuntos;k++){
				//~ cerr << indy[i][j][k] << " ";
			//~ }
			//~ cerr << endl;
		//~ }
	//~ }
//Fin de lectura de la instacia
///////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	
	
	//Se inicializa el ambiente de cplex
	
	env = CPXopenCPLEX(&status);
	
	if(env == NULL){
		char errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "No se pudo inicializar el ambiente.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		goto TERMINATE;
	}
	
	//Para habilitar el indicador por pantalla
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status){
		fprintf(stderr,"Fallo el indicador por pantalla, error %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
	if(status){
		fprintf(stderr,"No se pudo habilitar data checking, error %d\n", status);
		goto TERMINATE;
	}

	
	
	//Ahora se crea el problema
	
	lp = CPXcreateprob(env, &status, "formInicial");
	if(lp==NULL){
		fprintf(stderr, "No se pudo crear el problema, error %d\n", status);
		goto TERMINATE;
	}
	
	
	
	
	
	/* Turn on traditional search for use with control callbacks */

   status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
   if ( status )  goto TERMINATE;
   
   status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
   if(status) goto TERMINATE;

   /* Let MIP callbacks work on the original model */

   //~ status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
   //~ if ( status )  goto TERMINATE;

   /* Create user cuts for noswot problem */

   //~ status = makeusercuts (env, lp, &usercutinfo);
   //~ if ( status )  goto TERMINATE;

   /* Set up to use MIP usercut callback */
	
	
	

	
	
	
	
	
/////////////////////////////////////////////////////////////////////////
//Aca se llena el problema
	#define NUMCOLS cantPuntos*cantPuntos+cantPuntos*cantPuntos*cantPuntos
	#define NUMROWS 1
	#define NUMNZ cantPuntos*cantPuntos+cantPuntos*cantPuntos*cantPuntos
	status=0;
	
	double obj[NUMCOLS];
	double lb[NUMCOLS];
	double ub[NUMCOLS];
	char *colname[NUMCOLS];
	char coltype[NUMCOLS];
	int rmatbeg[NUMROWS];
	int rmatind[NUMNZ];
	double rmatval[NUMNZ];
	double rhs[NUMROWS];
	char sense[NUMROWS];
	char *rowname[NUMROWS];
	
	//Se setea el sentido de la funcion objetivo
	CPXchgobjsen(env, lp, CPX_MIN);
	
	//Se define la funcion objetivo, los bounds y los nombres de columnas
	for(int i=0;i<cantPuntos;i++){
		for(int j=0;j<cantPuntos;j++){
			obj[indx[i][j]] = alfa*dist[i][j];
			lb[indx[i][j]] = 0.0;
			ub[indx[i][j]] = 1.0;
			coltype[indx[i][j]]='B';
			for(int k=0;k<cantPuntos;k++){
				obj[indy[i][j][k]] = beta*ang[i][j][k];
				lb[indy[i][j][k]] = 0.0;
				ub[indy[i][j][k]] = 1.0;
				coltype[indy[i][j][k]]='B';
			}
		}
	}
	//Ahora esta agregando sin nombre a las columnas, la idea es ponerlo bien
	//~ status = CPXnewcols(env, lp, NUMCOLS, obj, lb, ub, coltype, colname);
	status = CPXnewcols(env, lp, NUMCOLS, obj, lb, ub, coltype, NULL);
	if(status)goto ENDFILL;
	
	
	//Agrega las restricciones para que de cada nodo solo se pueda salir una vez(y no a si mismo)
	for(int i=0;i<cantPuntos;i++){
		int nozero = cantPuntos-1;
		rmatbeg[0]=0;
		rhs[0] = 1.0;
		sense[0] = 'E';
		int count = 0;
		for(int j=0;j<cantPuntos;j++){
			if(j==i)continue;
			rmatval[count]=1.0;
			rmatind[count++]=indx[i][j];
		}
		//rowname[0] = 
		status = CPXaddrows(env, lp, 0, NUMROWS, nozero, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if(status)goto ENDFILL;
	}
	
	//Agrega las restricciones para que a cada nodo solo se pueda entrar una vez(y no desde si mismo)
	for(int j=0;j<cantPuntos;j++){
		int nozero = cantPuntos-1;
		rmatbeg[0]=0;
		rhs[0]=1.0;
		sense[0]='E';
		int count=0;
		for(int i=0;i<cantPuntos;i++){
			if(i==j)continue;
			rmatval[count]=1.0;
			rmatind[count++]=indx[i][j];
		}
		status = CPXaddrows(env, lp, 0, NUMROWS, nozero, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if(status)goto ENDFILL;
	}
	//~ cerr << "Cantidad de columnas: " << CPXgetnumcols(env,lp) << endl;
	//Agrega las restricciones para que un angulo sea elegido obligatoriamente si se eligieron las dos aristas
	for(int i=0;i<cantPuntos;i++){
		for(int j=0;j<cantPuntos;j++){
			if(i==j)continue;
			for(int k=0;k<cantPuntos;k++){
				//~ ////~ if(k==i||k==j)continue;
				if(k==j)continue;
				int nozero = 3;
				rmatbeg[0]=0;
				rhs[0]=1.0;
				sense[0]='L';
				rmatind[0]=indx[i][j];
				rmatval[0]=1.0;
				rmatind[1]=indx[j][k];
				rmatval[1]=1.0;
				rmatind[2]=indy[i][j][k];
				rmatval[2]=-1.0;
				//~ cerr << "Indice de la variable de angulo: " << rmatval[2] << endl;
				status = CPXaddrows(env, lp, 0, NUMROWS, nozero, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if(status)goto ENDFILL;
			}
		}
	}
	
ENDFILL:
	if(status){
		fprintf(stderr, "No se pudo llenar el problema, error %d\n", status);
		goto TERMINATE;
	}
	
/////////////////////////////////////////////////////////////////////////	
	//Optimizacion del problema
	
	
	//Aca quiero setear para hacer callbacks
	
	//cut callback con mi funcion
	status = CPXsetusercutcallbackfunc(env, mycutcallback, NULL);
	//cut cullback sin ninguna funcion
	//~ status = CPXsetusercutcallbackfunc(env, NULL, NULL);
	//~ status = CPXsetcutcallbackfunc(env, mycutcallback, NULL);
	if(status){
		cerr << "pincho callback" << endl;
		goto TERMINATE;
	}  
	
	
	
	
	
	status = CPXmipopt(env,lp);
	if(status){
		cout << "El status es: " << status << endl;
		fprintf(stderr,"Fallo la optimizacion\n");
		goto TERMINATE;
	}
/////////////////////////////////////////////////////////////////////////	
	//Pido la solucion y la muestro
	
	//Prueba para ver el status del problema
	//~ cerr << "STATUSSS: " << CPXgetstat(env,lp) << endl;
	
	//Para saber el valor de la funcion objetivo
	status = CPXgetobjval (env, lp, &objval);
	if ( status ) {
		fprintf (stderr,"No MIP objective value available.  Exiting...\n");
		goto TERMINATE;
	}
	
	printf ("Valor de la solucion  = %f\n\n", objval);
	
	
	
	
	
	cur_numrows=CPXgetnumrows(env,lp);
	cur_numcols=CPXgetnumcols(env,lp);
	
	x = (double *) malloc(cur_numcols * sizeof(double));
	slack = (double *) malloc(cur_numrows * sizeof(double));
	dj = (double *) malloc(cur_numcols * sizeof(double));
	pi = (double *) malloc(cur_numrows * sizeof(double));
	
	if(x==NULL||slack==NULL||dj==NULL||pi==NULL){
		status = CPXERR_NO_MEMORY;
		fprintf(stderr, "No se pudieron alocar las variables\n");
		goto TERMINATE;
	}
	
	
	status = CPXgetx (env, lp, x, 0, cur_numcols-1);
	if ( status ) {
		fprintf (stderr, "Failed to get optimal integer x.\n");
		goto TERMINATE;
	}

	status = CPXgetslack (env, lp, slack, 0, cur_numrows-1);
	if ( status ) {
		fprintf (stderr, "Failed to get optimal slack values.\n");
		goto TERMINATE;
	}

	//~ for(i=0;i<cur_numrows;i++){
		//~ printf("Row %d:  Slack = %10f\n",i,slack[i]);
	//~ }

	for(j=0;j<cur_numcols;j++) {
		printf("Column %d:  Value = %10f ",j,x[j]);
		cout << revind[j] << endl;
	}

   /* Ahora escribo el problema a un archivo */

   
	//~ 
	status = CPXwriteprob (env, lp, "mipex1.lp", NULL);
	if ( status ) {
		fprintf (stderr, "Failed to write LP to disk.\n");
		goto TERMINATE;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
///////////////////////////////////////////////////////////////////////////
// Todo este pedazo servia para tener informacion de la solucion, pero no sirve para MIP
	//Esta linea no anda para mipopt
	//~ status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
	//~ if(status){
		//~ fprintf(stderr, "No se pudo obtener la solucion\n");
		//~ goto TERMINATE;
	//~ }
	//~ 
	//~ 
	//~ printf("\nSolution status = %d\n", solstat);
	//~ printf("Valor de la solucion = %f\n\n", objval);
	//~ 
	//~ for(i=0;i<cur_numrows;i++){
		//~ printf("Fila %d: Slack = %10f Pi = %10f\n", i, slack[i], pi[i]);
	//~ }
	//~ for(j=0;j<cur_numcols;j++){
		//~ printf("Columna %d: Valor = %10f Costo reducido = %10f\n", j, x[j], dj[j]);
	//~ }
	
/////////////////////////////////////////////////////////////////////////
TERMINATE:
	puts("Termino");
	return status;
}




///////////////////////////////////////


///////////////////////////////////////


///////////////////////////////////////
static void free_and_null(char **ptr){
	if(*ptr!=NULL){
		free(*ptr);
		*ptr=NULL;
	}
}
///////////////////////////////////////
static int CPXPUBLIC 
mycutcallback (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               int        *useraction_p)
{
   cout << "Aca estoy haciendo callbacks!" << endl;
   //~ *useraction_p = CPX_CALLBACK_DEFAULT;
   return 0;
} /* END mycutcallback */


