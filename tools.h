//
// Created by josse on 7/15/2020.
//
#include <fstream>
#include <string>


void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
            case INT_FLOAT:
                int e0; float condition_value;
                file >> e0 >> condition_value;
                item_list[i].setValues(0,0,0,0,0,e0,condition_value,0);
                break;
            case INT_FLOAT_FLOAT_FLOAT:
                int e; float coord1, coord2, coord3;
                file >> e >> coord1 >> coord2 >> coord3;
                item_list[i].setValues(e,coord1,coord2,coord3,0,0,0,0);
                break;
            case INT_INT_INT_INT_INT:
                int element,node1,node2,node3,node4;
                file >> element >> node1 >> node2 >> node3 >> node4;
                item_list[i].setValues(element,0,0,0,node1,node2,node3,node4);
                break;
        }
    }
}

void correctConditions(int n,condition *list,int *indices){
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1();

    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}

void addExtension(char *newfilename,char *filename, string extension){
    int ori_length = strlen(filename);
    int ext_length = extension.length();
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void correctIndices(int n,condition* list,int total){
    for(int i=0;i<n;i++)
        list[i].setNode1(list[i].getNode1()+total);
}

void addArray(int* index, int n,condition* list, condition* listf){
    for(int i=0;i<n;i++){
        listf[*index] = list[i];
        (*index)++;
    }
}

void fusionCond(int n1,condition* list1,int n2,condition* list2,int n3,condition* list3,int n4,condition* list4,condition* list){
    int index = 0;
    addArray(&index,n1,list1,list);
    addArray(&index,n2,list2,list);
    addArray(&index,n3,list3,list);
    addArray(&index,n4,list4,list);
}

void fusionNeumann(int n1,condition* list1,int n2,condition* list2,int n3,condition* list3,int n4,condition* list4,condition* list){
    int index = 0;
    addArray(&index,n1,list1,list);
    addArray(&index,n2,list2,list);
    addArray(&index,n3,list3,list);
    addArray(&index,n4,list4,list);
}


void leerMallaYCondiciones(mesh &m,char *filename){
    char inputfilename[150];
    ifstream file;
    float kappa, cconst, omega, alpha, m_x, m_y, m_z, phi;
    int nnodes, neltos, ndirich_o, ndirich_m, ndirich_g, ndirich_t, nneu_o, nneu_m, nneu_g, nneu_t;
    condition *dirichlet_o, *dirichlet_m, *dirichlet_g, *dirichlet_t;
    condition *neumann_o, *neumann_m, *neumann_g, *neumann_t;

    addExtension(inputfilename,filename,".dat");
    file.open(inputfilename);

    cout << "leer constantes" << endl;
    file >> kappa >> cconst >> omega >> alpha >> m_x >> m_y >> m_z >> phi;

    cout << "leer numero de nodos, elementos, dirichlet" << endl;
    file >> nnodes >> neltos >> ndirich_o >> ndirich_m >> ndirich_g >> ndirich_t >> nneu_o >> nneu_m >> nneu_g >> nneu_t;

    m.setParameters(kappa, cconst, omega, alpha, m_x, m_y, m_z, phi);
    m.setSizes(nnodes,neltos,ndirich_o+ndirich_m+ndirich_g+ndirich_t, nneu_o + nneu_m + nneu_g + nneu_t);
    m.createData();

    dirichlet_o = new condition[ndirich_o];
    dirichlet_m = new condition[ndirich_m];
    dirichlet_g = new condition[ndirich_g];
    dirichlet_t = new condition[ndirich_t];

    //NEUMANN
    neumann_o = new condition[nneu_o];
    neumann_m = new condition[nneu_m];
    neumann_g = new condition[nneu_g];
    neumann_t = new condition[nneu_t];

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT_INT,m.getElements());


    obtenerDatos(file,DOUBLELINE,ndirich_o,INT_FLOAT,dirichlet_o);
    obtenerDatos(file,DOUBLELINE,ndirich_m,INT_FLOAT,dirichlet_m);
    obtenerDatos(file,DOUBLELINE,ndirich_g ,INT_FLOAT,dirichlet_g);
    obtenerDatos(file,DOUBLELINE,ndirich_t,INT_FLOAT,dirichlet_t);

    obtenerDatos(file,DOUBLELINE,nneu_o,INT_FLOAT,neumann_o);
    obtenerDatos(file,DOUBLELINE,nneu_m,INT_FLOAT,neumann_m);
    obtenerDatos(file,DOUBLELINE,nneu_g ,INT_FLOAT,neumann_g);
    obtenerDatos(file,DOUBLELINE,nneu_t,INT_FLOAT,neumann_t);

    file.close();

    //La primera condición no es necesario ajustarla
    correctIndices(ndirich_m,dirichlet_m, nnodes);
    correctIndices(ndirich_g,dirichlet_g,2*nnodes);
    correctIndices(ndirich_t,dirichlet_t,3*nnodes);

    //Corrigiendo indices de Neumann
    correctIndices(nneu_m,neumann_m, nnodes);
    correctIndices(nneu_g,neumann_g,2*nnodes);
    correctIndices(nneu_t,neumann_t,3*nnodes);


    fusionCond(ndirich_o,dirichlet_o,ndirich_m,dirichlet_m,ndirich_g,dirichlet_g,ndirich_t,dirichlet_t,m.getDirichlet());
    fusionCond(nneu_o,neumann_o,nneu_m,neumann_m,nneu_g,neumann_g,nneu_t,neumann_t,m.getNeumann());

    correctConditions(ndirich_o+ndirich_m+ndirich_g+ndirich_t,m.getDirichlet(),m.getDirichletIndices());
    correctConditions(nneu_o+nneu_m+nneu_g+nneu_t,m.getNeumann(),m.getNeumannIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

int getIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return i;
    return -1;
}

int getIndex(int v, int s, Vector vec){
    for(int i=0;i<s;i++)
        if(vec.at(i)==v) return i;
    return -1;
}

int *createNonDirichletIndices(int nn,int nd,int *dirich_indices){
    int *ndi = new int[4*nn-nd];
    int pos = 0;
    for(int i=1;i<=4*nn;i++)
        if(!findIndex(i,nd,dirich_indices)){
            ndi[pos] = i;
            pos++;
        }
    return ndi;
}

void writeResults(mesh m,Vector T,char *filename){
    char outputfilename[150];

    int nn = m.getSize(NODES);
    int nd = m.getSize(DIRICHLET);
    int *dirich_indices = m.getDirichletIndices();
    int *non_dirich_indices = createNonDirichletIndices(nn,nd,dirich_indices);
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";

    file << "Result \"Vector R\" \"Load Case 1\" 1 Vector OnNodes\nComponentNames \"o\" \"m\" \"g\"\nValues\n";

    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1,nd,dirich_indices);
        //Imprimiendo O
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << " ";
        else{
            int T_index = getIndex(i+1,4*nn-nd,non_dirich_indices);
            file << i+1 << " " << T.at(T_index) << " ";
        }
        //Imprimiendo M
        d_index = getIndex(i+1+nn,nd,dirich_indices);
        if(d_index != -1)
            file << dirich[d_index].getValue() << " ";
        else{
            int T_index = getIndex(i+1+nn,4*nn-nd,non_dirich_indices);
            file << T.at(T_index) << " ";
        }
        //Imprimiendo G
        d_index = getIndex(i+1+2*nn,nd,dirich_indices);
        if(d_index != -1)
            file << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+2*nn,4*nn-nd,non_dirich_indices);
            file << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file << "\nResult \"Escalar t\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"t\"\nValues\n";

    cout <<"Write results T" << endl;
    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1+3*nn,nd,dirich_indices);
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+3*nn,4*nn-nd,non_dirich_indices);

            file << i+1 << " " << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file.close();
}


