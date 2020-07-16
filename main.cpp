#include <iostream>
#include "classes.h"
#include "math_tools.h"
#include "tools.h"
#include "display_tools.h"
#include "sel_tools.h"
#include "components.h"
#include "sel.h"
#include "assembly.h"


int main(int argc, char *argv[])
{
    char filename[150];
    if(argv[1] != NULL){
        strcpy(filename,argv[1]);
    }else{
        char filename2[150] =  "Example";
        strcpy(filename,filename2);
    }


    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACION"<<"N DEL ME"<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- MODELO DE PROYECTO DE TSC01/2020\n" << "\t- 3 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";
    mesh m;
    cout << "leer" << endl;

    leerMallaYCondiciones(m,filename);
    cout << "terminar de leer" << endl;

    cout << "crear sistemas locales" << endl;
    crearSistemasLocales(m,localKs,localbs);

    cout << "ASM" << endl;
    zeroes(K,4*m.getSize(NODES));
    zeroes(b,4*m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);

    applyNeumann(m, b);
    applyDirichlet(m,K,b);

    cout << "K Global: " << endl;
    showMatrix(K);
    cout << endl;

    cout << "b Global: " << endl;
    showVector(b);
    cout << endl;

    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    writeResults(m,T,filename);

    return 0;
}

