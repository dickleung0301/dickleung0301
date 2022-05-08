#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;
// setting of the physical conditions of the problem
double node_list[3][2] = {{0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}};
int element_list[3][2] = {{1, 2}, {2, 3}, {3, 1}};
int boundary_condition[3][2] = {{-1, -1}, {1, -1}, {1, 1}};
double force[3][2] = { {0.0, 0.0}, {0.0, 0.0}, {0.0, -20.0}};
double displacement[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
double young_modulus = pow(10.0, 6);
double cross_sectional_area = 0.01;
// problem size
int problem_dim = sizeof(node_list[0])/sizeof(double);
int num_node = sizeof(node_list)/sizeof(node_list[0]);
int num_element = sizeof(element_list)/sizeof(element_list[0]);
int node_per_element = sizeof(element_list[0])/sizeof(int);
// define the using functions
void assign_BCs(double* extended_node_list, int &DOFs, int &DOCs){
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < problem_dim; j++){
            if(extended_node_list[i * 6 * problem_dim + j + problem_dim] == -1){
                DOCs -= 1;
                extended_node_list[i * 6 * problem_dim + j + 2 * problem_dim] = DOCs;
            }
            else{
                DOFs += 1;
                extended_node_list[i * 6 * problem_dim + j + 2 * problem_dim] = DOFs;
            }
        }
    }
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < problem_dim; j++){
            if(extended_node_list[i * 6 * problem_dim + j + problem_dim] < 0){
                extended_node_list[i * 6 * problem_dim + j + 3 * problem_dim] = abs(extended_node_list[i * 6 * problem_dim + j + 2 * problem_dim])  + DOFs;
            }
            else{
                extended_node_list[i * 6 * problem_dim + j + 3 * problem_dim] = abs(extended_node_list[i * 6 * problem_dim + j + 2 * problem_dim]);
            }
        }
    }
    DOCs = abs(DOCs);
}

void element_stiffness(int first_element, int second_element, double* element_siffness, double* extended_node_list, double young_modulus, double cross_sectional_area){
    // extract the coordinate of nodes
    double x1 = extended_node_list[(first_element - 1) * 6 * problem_dim];
    double y1 = extended_node_list[(first_element - 1) * 6 * problem_dim + 1];
    double x2 = extended_node_list[(second_element - 1) * 6 * problem_dim];
    double y2 = extended_node_list[(second_element - 1) * 6 * problem_dim + 1];
    // calculate the length, cosine and sine of the element
    double temp = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    double length = sqrt(temp);
    double c = (x2 - x1)/length;
    double s = (y2 - y1)/length;
    // calculate the stiffness of the element
    element_siffness[0] = (young_modulus * cross_sectional_area/length) * c * c;
    element_siffness[1] = (young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[2] = -(young_modulus * cross_sectional_area/length) * c * c;
    element_siffness[3] = -(young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[4] = (young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[5] = (young_modulus * cross_sectional_area/length) * s * s;
    element_siffness[6] = -(young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[7] = -(young_modulus * cross_sectional_area/length) * s * s;
    element_siffness[8] = -(young_modulus * cross_sectional_area/length) * c * c;
    element_siffness[9] = -(young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[10] = (young_modulus * cross_sectional_area/length) * c * c;
    element_siffness[11] = (young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[12] = -(young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[13] = -(young_modulus * cross_sectional_area/length) * s * s;
    element_siffness[14] = (young_modulus * cross_sectional_area/length) * c * s;
    element_siffness[15] = (young_modulus * cross_sectional_area/length) * s * s;
    cout << "Element stiffness : " << endl;
    for(int i = 0; i < 16; i ++)
    cout << element_siffness[i] << ' ';
    cout << endl;
}

void assemble_stiffness(double* stiffness_matrix, double* element_siffness, double* extended_node_list, double young_modulus, double cross_sectional_area){
    // define first element and second element
    int first_element = 0, second_element = 0;
    // calculate the element stiffness and assemble it to the global stiffness
    for(int i = 0; i < num_element; i++){
        first_element = element_list[i][0];
        second_element = element_list[i][1];
        element_stiffness(first_element, second_element, element_siffness, extended_node_list, young_modulus, cross_sectional_area);
        for(int j = 0; j < node_per_element; j++)
            for(int k = 0; k < problem_dim; k++)
                for(int l = 0; l < node_per_element; l++)
                    for(int m = 0; m < problem_dim; m++){
                        int temp, temp2;
                        if(j == 0) temp = first_element;
                        if(j == 1) temp = second_element;
                        if(l == 0) temp2 = first_element;
                        if(l == 1) temp2 = second_element;
                        int row = extended_node_list[(temp - 1) * 6 * problem_dim + k + 3 * problem_dim];
                        int column = extended_node_list[(temp2 - 1) * 6 * problem_dim + m + 3 * problem_dim];
                        double value = element_siffness[(j * problem_dim + k) * node_per_element * problem_dim + l * problem_dim + m];
                        stiffness_matrix[(row - 1) * num_node * problem_dim + column - 1] = stiffness_matrix[(row - 1) * num_node * problem_dim + column - 1] + value;
                    }
    }
}

int main()
{
    int DOFs = 0;
    int DOCs = 0;
    cout << "young's modulus: " << young_modulus << endl;
    // initialization of extended node list
    double* extended_node_list = new double[num_node * 6 * problem_dim];
    for(int i = 0; i < num_node; i++)
        for(int j = 0; j < 6 * problem_dim; j++)
            extended_node_list[i * 6 * problem_dim + j] = 0.0;
    // initialization of the node coordinates
    for(int i = 0; i < num_node; i++)
        for(int j = 0; j < problem_dim; j++)
            extended_node_list[i * 6 * problem_dim + j] = node_list[i][j];
    // initialization of whether the node can move or not
    for(int i = 0; i < num_node; i++)
        for(int j = problem_dim; j < 2 * problem_dim; j++)
            extended_node_list[i * 6 * problem_dim + j] = boundary_condition[i][j - problem_dim];
    //assigning the boundary conditions
    assign_BCs(extended_node_list, DOFs, DOCs);
    // checking extended node list
    cout << "The original version of extended node list:" << endl;
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < 6 * problem_dim; j++) 
            cout << extended_node_list[i * 6 * problem_dim + j] << ",";
        cout << endl;
    }
    // initialization of stiffness matrix and element stiffness matrix
    double* stiffness_matrix = new double[num_node * problem_dim * num_node * problem_dim];
    double* element_siffness = new double[node_per_element * problem_dim * node_per_element * problem_dim];
    // calculate the global stiffness
    assemble_stiffness(stiffness_matrix, element_siffness, extended_node_list, young_modulus, cross_sectional_area);
    // checking global stiffness matrix
    cout << "The stiffness matrix:" << endl;
    for(int i = 0; i < num_node * problem_dim; i++){
        for(int j = 0; j < num_node * problem_dim; j++) 
            cout << stiffness_matrix[i * num_node * problem_dim + j] << ",";
        cout << endl;
    }
    // initialization of force and displacement
    for(int i = 0; i < num_node; i++)
        for(int j = 0; j < problem_dim; j++)
            extended_node_list[i * 6 * problem_dim + 4 * problem_dim + j] = displacement[i][j];
    for(int i = 0; i < num_node; i++)
        for(int j = 0; j < problem_dim; j++)
            extended_node_list[i * 6 * problem_dim + 5 * problem_dim + j] = force[i][j];
    // checking extended node list after initialization of force and displacement
    cout << "The updated version of extended node list:" << endl;
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < 6 * problem_dim; j++) 
            cout << extended_node_list[i * 6 * problem_dim + j] << ",";
        cout << endl;
    }
    // assemble force
    vector<double> forcep;
    int DOF = 0;
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < problem_dim; j++){
            if(extended_node_list[i * 6 * problem_dim + j + problem_dim] == 1){
                DOF += 1;
                forcep.push_back(extended_node_list[i * 6 * problem_dim + j + 5 * problem_dim]);
            }
        }
    }
    // assemble displacement 
    vector<double> displacementp;
    int DOC = 0;
    for(int i = 0; i < num_node; i++){
        for(int j = 0; j < problem_dim; j++){
            if(extended_node_list[i * 6 * problem_dim + j + problem_dim] == 1){
                DOC += 1;
                displacementp.push_back(extended_node_list[i * 6 * problem_dim + j + 4 * problem_dim]);
            }
        }
    }
    // calculate the unknown force and displacement
    cout << "K_UU :" << endl;
    double* K_UU = new double[DOFs * DOFs];
    for(int i = 0; i < DOFs; i++){
        for(int j = 0; j < DOFs; j++){
            K_UU[i * DOFs + j] = stiffness_matrix[i * num_node * problem_dim + j];
            cout << K_UU[i * DOFs + j] << ' ';
        }
        cout << endl;
    }
    cout << "K_UP :" << endl;
    double* K_UP = new double[DOFs * DOCs];
    for(int i = 0; i < DOFs; i++){
        for(int j = 0; j < DOCs; j++){
            K_UP[i * DOCs + j] = stiffness_matrix[i * num_node * problem_dim + DOFs + j];
            cout << K_UP[i * DOCs + j] << ' ';
        }
        cout << endl;
    }
    cout << "K_PU :" << endl;
    double* K_PU = new double[DOFs * DOCs];
    for(int i = 0; i < DOCs; i++){
        for(int j = 0; j < DOFs; j++){
            K_PU[i * DOFs + j] = stiffness_matrix[(i + DOFs) * num_node * problem_dim + j]; 
            cout << K_PU[i * DOFs + j] << ' ';
        }
        cout << endl;
    }
    cout << "K_PP :" << endl; 
    double* K_PP = new double[DOFs * DOCs];
    for(int i = 0; i < DOCs; i++){
        for(int j = 0; j < DOFs; j++){
            K_PP[i * DOFs + j] = stiffness_matrix[(i + DOFs) * num_node * problem_dim + DOFs + j];
            cout << K_PP[i * DOFs + j] << ' ';
        }
        cout << endl;
    }
    double determinant = K_UU[0] * (K_UU[4] * K_UU[8] - K_UU[5] * K_UU[7]) 
                        - K_UU[1] * (K_UU[3] * K_UU[8] - K_UU[5] * K_UU[6]) 
                        + K_UU[2] * (K_UU[3] * K_UU[7] - K_UU[4] * K_UU[6]);
    cout << "determinant :" << determinant << endl;
    cout << "inverse_K_UU" << endl;
    double* inverse_K_UU = new double[DOFs * DOFs];
    inverse_K_UU[0] = (K_UU[4] * K_UU[8] - K_UU[5] * K_UU[7]) / determinant;
    inverse_K_UU[1] = (K_UU[2] * K_UU[7] - K_UU[1] * K_UU[8]) / determinant;
    inverse_K_UU[2] = (K_UU[1] * K_UU[5] - K_UU[2] * K_UU[4]) / determinant;
    inverse_K_UU[3] = (K_UU[5] * K_UU[6] - K_UU[3] * K_UU[8]) / determinant;
    inverse_K_UU[4] = (K_UU[0] * K_UU[8] - K_UU[2] * K_UU[6]) / determinant;
    inverse_K_UU[5] = (K_UU[2] * K_UU[3] - K_UU[0] * K_UU[5]) / determinant;
    inverse_K_UU[6] = (K_UU[3] * K_UU[7] - K_UU[4] * K_UU[6]) / determinant;
    inverse_K_UU[7] = (K_UU[1] * K_UU[6] - K_UU[0] * K_UU[7]) / determinant;
    inverse_K_UU[8] = (K_UU[0] * K_UU[4] - K_UU[1] * K_UU[3]) / determinant;
    for(int i = 0; i < DOFs; i++){
        for(int j = 0; j < DOFs; j++){
            cout << inverse_K_UU[i * DOFs + j] << ' ';
        }
        cout << endl;
    }
    cout << "size of forcep: " << forcep.size() << endl;
    cout << "size of displacementp: " << displacementp.size() << endl;
    cout << "DOFs :" << DOFs << endl;
    cout << "DOCs :" << DOCs << endl;

    double* F = new double[DOFs];
    for(int i = 0; i < DOFs; i++){
        for(int j = 0; j < DOCs; j++)
            F[i] += K_UP[i * DOCs + j] * displacementp[j];
    }

    for(int i = 0; i < DOFs; i++){
        F[i] = forcep[i] - F[i];
    }

    double* U_u = new double[DOFs];
    for(int i = 0; i < DOFs; i++){
        for(int j = 0; j < DOFs; j++)
            U_u[i] += inverse_K_UU[i * DOFs + j] * F[j];
    }

    double* temp = new double[DOCs];
    for(int i = 0; i < DOCs; i++){
        for(int j = 0; j < DOFs; j++)
            temp[i] += K_PU[i * DOFs + j] * U_u[j];
    }

    double* temp2 = new double[DOCs];
    for(int i = 0; i < DOCs; i++){
        for(int j = 0; j < DOCs; j++)
            temp2[i] += K_PP[i * DOCs + j] * displacementp[j];
    }

    double* F_u = new double[DOCs];
    for(int i = 0; i < DOCs; i++)
        F_u[i] = temp[i] + temp2[i];
    
    cout << "The force: ";
    for(int i = 0; i < DOCs; i++)
        cout << F_u[i] << ' ';
    cout << endl;

    cout << "The displacement: ";
    for(int i = 0; i < DOFs; i++)
        cout << U_u[i] << ' ';
    cout << endl;

    // free the dynamically allocated memory
    free(extended_node_list);
    free(stiffness_matrix);
    free(element_siffness);
    free(K_UU);
    free(K_UP);
    free(K_PU);
    free(K_PP);
    free(inverse_K_UU);
    free(F);
    free(U_u);
    free(temp);
    free(temp2);
    free(F_u);
}