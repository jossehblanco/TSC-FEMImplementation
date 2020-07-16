//
// Created by josse on 7/15/2020.
//

enum lines {NOLINE,SINGLELINE,DOUBLELINE};
enum modes {NOMODE,INT_FLOAT,INT_FLOAT_FLOAT_FLOAT,INT_INT_INT_INT_INT};
enum parameters {KAPPA, CCONST, OMEGA, ALPHA, M_X, M_Y, M_Z,  PHI};
enum sizes {NODES,ELEMENTS,DIRICHLET, NEUMANN};
enum coords {EQUIS,YE,ZETA};

class item{
protected:
    int id;
    float x;
    float y;
    float z;
    int node1;
    int node2;
    int node3;
    int node4;
    float value;
public:
    void setId(int identifier) {
        id = identifier;
    }

    void setX(float x_coord) {
        x = x_coord;
    }

    void setY(float y_coord) {
        y = y_coord;
    }

    void setZ(float z_coord) {
        z = z_coord;
    }

    void setNode1(int node_1) {
        node1 = node_1;
    }

    void setNode2(int node_2) {
        node2 = node_2;
    }

    void setNode3(int node_3) {
        node3 = node_3;
    }

    void setNode4(int node_4) {
        node4 = node_4;
    }

    void setValue(float value_to_assign) {
        value = value_to_assign;
    }

    int getId() {
        return id;
    }

    float getX() {
        return x;
    }

    float getY() {
        return y;
    }

    float getZ() {
        return z;
    }

    int getNode1() {
        return node1;
    }

    int getNode2() {
        return node2;
    }

    int getNode3() {
        return node3;
    }

    int getNode4() {
        return node4;
    }

    float getValue() {
        return value;
    }

    virtual void setValues(int a,float b,float c,float d,int e,int f,float g, int h)=0;


};

class node: public item{

public:
    void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
        id = a;
        x = b;
        y = c;
        z = d;
    }

};

class element: public item{

public:
    void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
        id = a;
        node1 = e;
        node2 = f;
        node3 = g;
        node4 = h;
    }

};

class condition: public item{

public:

    void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
        node1 = f;
        value = g;
    }

};

class mesh{
    float parameters[6];
    int sizes[4];
    node *node_list;
    element *element_list;
    int *indices_dirich;
    int *indices_neumann;
    condition *dirichlet_list;
    condition *neumann_list;
public:
    //KAPPA, CCONST, OMEGA, ALPHA, M_X, M_Y, M_Z,  PHI
    void setParameters(float kappa,float cconst, float omega, float alpha, float m_x, float m_y, float m_z, float phi){
        parameters[KAPPA] = kappa;
        parameters[CCONST] = cconst;
        parameters[OMEGA] = omega;
        parameters[ALPHA] = alpha;
        parameters[M_X] = m_x;
        parameters[M_Y] = m_y;
        parameters[M_Z] = m_z;
        parameters[PHI] = phi;

    }
    void setSizes(int nnodes,int neltos,int ndirich, int nneu){
        sizes[NODES] = nnodes;
        sizes[ELEMENTS] = neltos;
        sizes[DIRICHLET] = ndirich;
        sizes[NEUMANN] = nneu;
    }
    int getSize(int s){
        return sizes[s];
    }
    float getParameter(int p){
        return parameters[p];
    }
    void createData(){
        node_list = new node[sizes[NODES]];
        element_list = new element[sizes[ELEMENTS]];
        indices_dirich = new int[sizes[DIRICHLET]];
        dirichlet_list = new condition[sizes[DIRICHLET]];
        indices_neumann = new int[sizes[NEUMANN]];
        neumann_list = new condition[sizes[NEUMANN]];
    }
    node* getNodes(){
        return node_list;
    }
    element* getElements(){
        return element_list;
    }
    int* getDirichletIndices(){
        return indices_dirich;
    }
    condition* getDirichlet(){
        return dirichlet_list;
    }
    int* getNeumannIndices(){
        return indices_neumann;
    }
    condition* getNeumann(){
        return neumann_list;
    }
    node getNode(int i){
        return node_list[i];
    }
    element getElement(int i){
        return element_list[i];
    }
    condition getCondition(int i, int type){
        switch(type){
            case DIRICHLET:
                return dirichlet_list[i];
            case NEUMANN:
                return neumann_list[i];
        }

    }
};


