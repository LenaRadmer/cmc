#if 0
#include <../example/common/cmc_t8_refine_coarsen_mesh.h>
#include <iostream>
#include <t8.h>                                 /*General t8code header. */
#include "t8_cmesh.h"                           /*cmesh definition and basic interface. */
#include <t8_cmesh_vtk_writer.h>                /*cmesh-writer interface. */
#include "t8_cmesh/t8_cmesh_examples.h"         /*collection of exemplary cmeshes. */
#include <t8_forest/t8_forest_general.h>        /* forest definition and general interface. */
#include <t8_forest/t8_forest_io.h>             /* forest io interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>    /*geometrical information of the forest*/
#include <t8_vec.h>                             /*Basic operations on 3D vectors*/
#include <t8_element_c_interface.h>
#include <t8_element.h>



/*Our search query, a particle together with a flag*/ 
typedef struct {
    double coordinates[3];              /*The coordinates of our particle*/
    int is_inside_partition;            /*Will be set to true if the particles lies inside this process' parallel partition*/

} t8_tutorial_search_particle_t;        /*Frage: Was ist das genau? Antwort: der Alias*/

/*Additional user data that we process during search*/
typedef struct {
    sc_array *particles_per_element;            /*For each element the number of particles inside it*/
    t8_locidx_t num_elements_searched;          /*The total number of elemennts created*/

} t8_tutorial_search_user_data_t;

/*The search callback. It will be called once per element and decides if we should continue with the children of the element or not.
 * We return 1, since we will continue as long as there are particles left to consider. It will stop when no queries are active or the element is a leaf element.
 * We increase a counter by one to count how many elements we looked at.*/
static int
t8_search_callback (t8_forest_t forest,
                    t8_locidx_t ltreeid,
                    const t8_element_t *element,
                    const int is_leaf, 
                    t8_element_array_t *leaf_elements,
                    t8_locidx_t tree_leaf_index,
                    void *query, 
                    size_t query_index) {
    T8_ASSERT (query == NULL);

    /*Get a pointer to our user data and increase the counter of searched elements*/
    t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
    T8_ASSERT (user_data != NULL);
    user_data->num_elements_searched++;
    return 1;
    }

/*The query callback. It will be called for each element once per active query. The return value determines wheter or not 
 *the query remains active for the children. If isLeaf = true, the element is a actual leaf element.
 *If isLeaf = true and the element contains a particle we will increase the counter for this element by one.
 *These counters are provided in an sc_arry as user data of the input forest.*/
static int
t8_search_query_callback (t8_forest forest,
                          t8_locidx_t ltreeid,
                          const t8_element_t *element,
                          const int is_leaf,
                          t8_element_array_t *leaf_elements,
                          t8_locidx_t tree_leaf_index,
                          void *query,
                          size_t query_index) {

    int particle_is_inside_element;
    /*Cast the query pointer to a particle pointer*/
    t8_search_particle_t *particle = (t8_search_particle_t *) query;
    }
#endif











#if 0
// Example program
#include <iostream>
#include <string>
#include <vector>
//navarrete@lrz.de

std::vector<int> matrix2vector(std::vector<std::vector<int>> matrix) {
    
    int rows = matrix.size();
    int colums = matrix[0].size();
    std::vector<int> serializedMatrix;
    for (int n=0; n<rows; n++) {
        for (int m=0; m<colums; m++) 
           serializedMatrix.push_back(matrix[n][m]);
    }
    return serializedMatrix;

}

std::vector<std::vector<int>> vector2matrix (std::vector<int> vec, int n, int m) {
    
    std::vector< std::vector<int> > matrix;
    for (int i=0; i<n; i++) {
        std::vector<int> zeile;
        for (int j=0; j<m; j++) 
            zeile.push_back(i*m+j);
        matrix.push_back(zeile);
    }
    return matrix;
    
}

//isOnes
/*
class Animal {
        int d;
        public:
        virtual int getD(){returns d} = 0;
        virtual void eat()=0; // abstract base class
        
        };
        
        class Dog: public Animal { public: int d; void eat() override; int getD()...};
        
        void Dog::eat() { std::cout << "eating meat" << std::endl; }
        
        class Cat: public Animal { public: void eat() override;};
        
        void Cat::eat() { std::cout << "eating fish" << std::endl; }
        
        class Cow: public Animal { public: void eat() override; };
        
        void Cow::eat() { std::cout << "eating grass" << std::endl; }
        
        
        std::vector<Animal *> animals;
        
        
        
        
        Animal *dog = new Dog();
        Animal *cat = new Cat();
        Animal *cow = new Cow();
        animals.push_back(dog);
        animals.push_back(cat);
        animals.push_back(cow);
        for(auto animalPtr: animals) {
            animalPtr->eat();
            }
            
        
*/


#include <iostream>
#include <vector>
int main () {
   /* std::cout << "Hello World!" << std::endl;
    int n = 1;          //declaration statement
    n = n + 1;          //expression
    std::cout << "n = " << n << '\n';


    return 0;
    */
    int N = 6;
    int M = 5;
  
    std::vector< std::vector<int> > matrix;
    for (int n=0; n<N; n++) {
        std::vector<int> zeile;
        for (int m=0; m<M; m++) 
            zeile.push_back(n*M+m);
        matrix.push_back(zeile);
    }
  
    std::vector<int> newMatrix;
    newMatrix = matrix2vector(matrix);
    std::cout <<"NewMatrix: ";
    for (int i=0; i<N*M; i++) 
        std::cout <<newMatrix[i] << " ";
    std::cout << std::endl;
  
    std::vector<std::vector<int>> matrixFromVector;
    matrixFromVector = vector2matrix(newMatrix, N, M);
    std::cout <<"matrixFromVector: ";
    for (int n=0; n<N; n++) {
        std::cout << " \n ";
            for (int m=0; m<M; m++) 
                std::cout <<matrixFromVector[n][m] << " ";
        }
    std::cout << std::endl;   
}
#endif














// Example program
#include <iostream>
#include <cmath>
#include <vector>
//egelhofer@lrz.de

using namespace std;

class Star {
    public: 
    const double mass;
    const double temperature;
    const double density;
    Star(double &m, double &t, double &rho) : mass(m), temperature(t), density(rho) {}
    virtual double lum() {double lum= pow(mass, 3.5) / temperature; return lum; };
    virtual ~Star() {};
    
};
class BlackHole: public Star { public: BlackHole(double &m, double &t, double &rho):Star(m, t, rho) {}; double lum() override; ~BlackHole() {}; };

double BlackHole::lum() {
    double lum = mass * pow(temperature, 4);
    return lum;
};

class WhiteDwarf: public Star { public: WhiteDwarf(double &m, double &t, double &rho):Star(m, t, rho) {}; double lum() override; ~WhiteDwarf() {}; };

double WhiteDwarf::lum() {
    double lum = density * mass;
    return lum;
};


        
int main()
{   
    int it = 1;
    double mass = 1.2;
    double temperature = 4.2;
    double density = 12.;
    Star *newBlackHole = new BlackHole(mass, temperature, density);
    cout << newBlackHole->lum() << endl;
    Star *newWhiteDwarf = new WhiteDwarf(mass, temperature, density);
    cout << newWhiteDwarf->lum() << endl;
    Star *newStar = new Star(mass, temperature, density);
    cout << newStar->lum() << endl;
    vector<Star*> list_of_stellar_types;
    list_of_stellar_types.push_back(newBlackHole);
    list_of_stellar_types.push_back(newWhiteDwarf);
    list_of_stellar_types.push_back(newStar);

    for (auto i: list_of_stellar_types) {
        cout << "Luminosity of  " << it <<" : " << i->lum() << endl;
        ++it;

    }
    return 0;
}
