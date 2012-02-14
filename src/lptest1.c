#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>

static int populatebyrow (CPXENVptr env, CPXLPptr lp);
static void free_and_null (char** ptr);

int main(int argc, char **argv){
	
	int solstat;
	double objval;
	double *x = NULL;
	double *pi = NULL;
	double *slack = NULL;
	double *dj = NULL;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status = 0;
	int i, j;
	int cur_numrows, cur_numcols;
	
	if(argc != 1){
		puts("Se usa sin argumentos");
		goto TERMINATE;
	}
	
	//Se inicializa el ambiente de CPLEX
	
	env = CPXopenCPLEX (&status);
	
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
	
	//Para habilitar data checking
	status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
	if(status){
		fprintf(stderr, "No se pudo habilitar data checking, error %d\n", status);
		goto TERMINATE;
	}
	
	
	//Crear el problema
	
	lp = CPXcreateprob(env, &status, "lptest1");
	if(lp == NULL){
		fprintf(stderr, "No se pudo crear el problema, error %d\n", status);
		goto TERMINATE;
	}
	
	//Ahora se llena el problema
	status = populatebyrow(env, lp);
	if(status){
		fprintf(stderr, "No se pudo llenar el problema, error %d\n", status);
		goto TERMINATE;
	}
	
	//Optimizar el problema
	status = CPXlpopt(env,lp);
	if(status){
		fprintf(stderr, "Fallo la optimizacion\n");
		goto TERMINATE;
	}
	
	//En vez de hardcodear los tamaños, se los preguntamos al problema
	cur_numrows = CPXgetnumrows(env,lp);
	cur_numcols = CPXgetnumcols(env,lp);
	
	x = (double *) malloc(cur_numcols * sizeof(double));
	slack = (double *) malloc(cur_numrows * sizeof(double));
	dj = (double *) malloc(cur_numcols * sizeof(double));
	pi = (double *) malloc(cur_numrows * sizeof(double));
	
	if(x==NULL||slack==NULL||dj==NULL||pi==NULL){
		status = CPXERR_NO_MEMORY;
		fprintf(stderr, "No se pudieron alocar las variables\n");
		goto TERMINATE;
	}

	status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
	if(status){
		fprintf(stderr, "No se pudo obtener la solucion\n");
		goto TERMINATE;
	}
	
	
	//Escribimos la solucion por pantalla
	
	printf("\nSolution status = %d\n", solstat);
	printf("Valor de la solucion = %f\n\n", objval);
	
	for(i=0;i<cur_numrows;i++){
		printf("Fila %d: Slack = %10f Pi = %10f\n", i, slack[i], pi[i]);
	}
	for(j=0;j<cur_numcols;j++){
		printf("Columna %d: Valor = %10f Costo reducido = %10f\n", j, x[j], dj[j]);
	}
	
	//Guardo el lp en un archivo
	status = CPXwriteprob(env,lp, "lptest1.lp", NULL);
	if(status){
		fprintf(stderr, "Fallo la escritura del lp a archivo\n");
		goto TERMINATE;
	}
	
TERMINATE:

	free_and_null((char**) &x);
	free_and_null((char**) &slack);
	free_and_null((char**) &pi);
	free_and_null((char**) &dj);
	
	if(lp!=NULL){
		status = CPXfreeprob(env,&lp);
		if(status){
			fprintf(stderr, "No se pudo liberar el problema, error %d\n", status);
		}
	}
	
	if(env != NULL){
		status = CPXcloseCPLEX(&env);
		if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "No se pudo cerrar el ambiente\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
		}
	}
	return status;
}


///////////////////////////////////////
static void free_and_null(char **ptr){
	if(*ptr!=NULL){
		free(*ptr);
		*ptr=NULL;
	}
}
///////////////////////////////////////


//Llenar por fila
static int populatebyrow(CPXENVptr env, CPXLPptr lp){
	#define NUMROWS 3
	#define NUMCOLS 3
	#define NUMNZ 9
	int		status = 0;
	double	obj[NUMCOLS];
	double	lb[NUMCOLS];
	double	ub[NUMCOLS];
	char  	*colname[NUMCOLS];
	int		rmatbeg[NUMROWS];
	int		rmatind[NUMNZ];
	double	rmatval[NUMNZ];
	double	rhs[NUMROWS];
	char 	sense[NUMROWS];
	char	*rowname[NUMROWS];
	
	//Primero cambio el sentido de la funcion objetivo
	CPXchgobjsen (env, lp, CPX_MAX);
	
	//Defino la funcion objetivo
	
	obj[0]=1.0;
	obj[1]=1.0;
	obj[2]=1.0;
	
	//Defino los bounds
	
	lb[0] = 0.0;
	lb[1] = 0.0;
	lb[2] = 0.0;
	ub[0] = CPX_INFBOUND;
	ub[1] = CPX_INFBOUND;
	ub[2] = CPX_INFBOUND;
	
	//Nombres de las columnas
	
	colname[0] = "x1";
	colname[1] = "x2";
	colname[2] = "x3";
	
	status = CPXnewcols(env, lp, NUMCOLS, obj, lb, ub, NULL, colname);
	if(status)goto TERMINATE;
	
	
	//Ahora se crean los constraints
	
	//Primera fila
	rmatbeg[0]=0;
	rowname[0]="c1";
	rmatind[0]=0;
	rmatval[0]=1.0;
	rmatind[1]=1;
	rmatval[1]=1.0;
	rmatind[2]=2;
	rmatval[2]=1.0;
	sense[0]='L';
	rhs[0]=20.0;
	
	//Segunda fila
	rmatbeg[1]=3;
	rowname[1]="c2";
	rmatind[3]=0;
	rmatval[3]=2.0;
	rmatind[4]=1;
	rmatval[4]=2.0;
	rmatind[5]=2;
	rmatval[5]=2.0;
	sense[1]='L';
	rhs[1]=40.0;
	
	//Segunda fila
	rmatbeg[2]=6;
	rowname[2]="c3";
	rmatind[6]=0;
	rmatval[6]=3.0;
	rmatind[7]=1;
	rmatval[7]=3.0;
	rmatind[8]=2;
	rmatval[8]=3.0;
	sense[2]='L';
	rhs[2]=60.0;
	
	status = CPXaddrows(env, lp, 0, NUMROWS, NUMNZ, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
	if(status)goto TERMINATE;
	
TERMINATE: 
	return status;
}
//////////////////////////////
