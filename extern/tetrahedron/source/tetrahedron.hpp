#include <cstring>

using namespace std;

namespace tetrahedron
{
    int i4_gcd ( int i, int j );
    int i4_lcm ( int i, int j );
    int i4_max ( int i1, int i2 );
    int i4_min ( int i1, int i2 );
    void i4vec_copy ( int n, int a1[], int a2[] );
    int i4vec_lcm ( int n, int v[] );
    double polygon_area_3d ( int n, double v[], double normal[] );
    double r8_acos ( double c );
    void r8_swap ( double *x, double *y );
    double r8_uniform_01 ( int &seed );
    void r8mat_copy ( int m, int n, double a1[], double a2[] );
    double r8mat_det_4d ( double a[] );
    int r8mat_solve ( int n, int rhs_num, double a[] );
    void r8mat_transpose_print ( int m, int n, double a[], string title );
    void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
      int ihi, int jhi, string title );
    double r8vec_angle_3d ( double u[], double v[] );
    void r8vec_copy ( int n, double a1[], double a2[] );
    double *r8vec_cross_3d ( double v1[3], double v2[3] );
    double r8vec_dot ( int n, double a1[], double a2[] );
    double r8vec_length ( int dim_num, double x[] );
    double r8vec_max ( int n, double r8vec[] );
    double r8vec_norm ( int n, double a[] );
    void r8vec_print ( int n, double a[], string title );
    void r8vec_transpose_print ( int n, double a[], string title );
    void r8vec_zero ( int n, double a[] );
    void shape_print ( int point_num, int face_num, int face_order_max,
      double point_coord[], int face_order[], int face_point[] );
    double *tetrahedron_barycentric ( double tetra[3*4], double p[3] );
    double *tetrahedron_centroid ( double tetra[3*4] );
    void tetrahedron_circumsphere ( double tetra[3*4], double &r, double pc[3] );
    bool tetrahedron_contains_point ( double tetra[3*4], double p[3] );
    double *tetrahedron_dihedral_angles ( double tetra[] );
    double *tetrahedron_edge_length ( double tetra[3*4] );
    void tetrahedron_edges ( double tetra[3*4], double ab[], double ac[],
      double ad[], double bc[], double bd[], double cd[] );
    void tetrahedron_face_angles ( double tetra[], double angles[] );
    void tetrahedron_face_areas ( double tetra[], double areas[] );
    void tetrahedron_insphere ( double tetra[3*4], double &r, double pc[3] );
    void tetrahedron_lattice_layer_point_next ( int c[], int v[], bool *more );
    void tetrahedron_lattice_point_next ( int c[], int v[], bool *more );
    double tetrahedron_quality1 ( double tetra[3*4] );
    double tetrahedron_quality2 ( double tetra[3*4] );
    double tetrahedron_quality3 ( double tetra[3*4] );
    double tetrahedron_quality4 ( double tetra[3*4] );
    void tetrahedron_rhombic_shape ( int point_num, int face_num,
      int face_order_max, double point_coord[], int face_order[],
      int face_point[] );
    void tetrahedron_rhombic_size ( int *point_num, int *edge_num,
      int *face_num, int *face_order_max );
    void tetrahedron_sample ( double tetra[3*4], int n, int &seed, double p[] );
    void tetrahedron_shape ( int point_num, int face_num, int face_order_max,
      double point_coord[], int face_order[], int face_point[] );
    void tetrahedron_size ( int *point_num, int *edge_num, int *face_num,
      int *face_order_max );
    double *tetrahedron_solid_angles ( double tetra[] );
    int tetrahedron_unit_lattice_point_num ( int s );
    double tetrahedron_volume ( double tetra[3*4] );
    void triangle_angles_3d ( double t[3*3], double angle[3] );
    double triangle_area_3d ( double t[3*3] );
}
