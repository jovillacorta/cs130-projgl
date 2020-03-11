#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color=0;
    state.image_depth=0;

    state.image_color = new pixel[state.image_height * state.image_width];

    for (unsigned int j = 0; j < state.image_height; j++){
    	for (unsigned int i = 0; i < state.image_width; i++){
    		state.image_color[i + state.image_width * j] = make_pixel(0,0,0);
    	}
    }

//    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_vertex a[3];
    data_geometry b[3];

	switch (type){
		case render_type::triangle: { // case for rendering triangle

            int num_Triangles = state.num_vertices / 3;

            for (int i = 0; i < num_Triangles; i++){
                for (int j = 0; j < 3; j++){
                    a[j].data = state.vertex_data + i * 3 * state.floats_per_vertex 
                                            + j * state.floats_per_vertex;
                    b[j].data= a[j].data;

                    state.vertex_shader(a[j], b[j], state.uniform_data);
                }

                const data_geometry *triangle[3] = {&b[0], &b[1], &b[2]};
                rasterize_triangle(state, triangle);
            }
            break;
        }
        
        case render_type::indexed: {
            break;
        }
        case render_type::fan: {
            break;
        }
        case render_type::strip: {
            break;
        }
			

		default: {
        }//
	}

//delete [] dataTriangle;
//    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{

    float i[3], j[3];

    //calculates i&j coordinates of vertices
    for (int q = 0; q < 3; ++q){
        float adjPosition;
        adjPosition = (in[q]->gl_Position[0]) / (in[q]->gl_Position[3]);
        i[q] = ((adjPosition + 1.0f) / 2.0f * state.image_width) - .5f;
        adjPosition = (in[q]->gl_Position[1]) / (in[q]->gl_Position[3]); 
        j[q] = ((adjPosition + 1.0f) / 2.0f * state.image_height) - .5f;       
    }

    // float min_x_Triangle = std::min(i[0],i[1],i[2]);
    // float min_y_Triangle = std::min(j[0],j[1],j[2]);
    // float max_x_Triangle = std::max(i[0],i[1],i[2]);
    // float max_y_Triangle = std::max(j[0],j[1],j[2]);

    // note i is i coordinate
    //      j is j coordinate
    //      index 0,1,2 is vertice A,B,C respectively

    const int A = 0;
    const int B = 1; 
    const int C = 2;

    float area_ABC = (i[A]*(j[B]-j[C]) + i[B]*(j[C]-j[A]) + i[C]*(j[A]-j[B])) * .5f;


    float alpha, beta, gamma;
    float k_alpha[3], k_beta[3], k_gamma[3];


    // calculate K values
    k_alpha[0] = (i[B] * j[C] - i[C] * j[B]) * .5f / area_ABC;
    k_alpha[1] = (j[B] - j[C]) * .5f / area_ABC;
    k_alpha[2] = (i[C] - i[B]) * .5f / area_ABC;

    k_beta[0] = (i[C] * j[A] - i[A] * j[C]) * .5f / area_ABC;
    k_beta[1] = (j[C] - j[A]) * .5f / area_ABC;
    k_beta[2] = (i[A] - i[C]) * .5f / area_ABC;

    k_gamma[0] = (i[A] * j[B] - i[B] * j[A]) * .5f / area_ABC;
    k_gamma[1] = (j[A] - j[B]) * .5f / area_ABC;
    k_gamma[2] = (i[B] - i[A]) * .5f / area_ABC;
    // pseudo code: for all x do 
    //               for all y do
    //                      compute alpha, beta, gamma for (x,y)
    //                          if (0<= alpha,beta,gamma<=1)
    //                              then c = alpha*ca + beta*cb + gamma*cc
    //                                 draw pixel (x,y) with color c
    for (unsigned int y = 0; y < state.image_height; y++){
        for (unsigned int x = 0; x < state.image_width; x++){
            
            for (int v=0; v < 3; v++){
                alpha = k_alpha[0] + k_alpha[1] * x + k_alpha[2] * y;
                beta = k_beta[0] + k_beta[1] * x + k_beta[2] * y;
                gamma = k_gamma[0] + k_gamma[1] * x + k_gamma[2] * y;
                // alpha = (0.5f * ((i[1] * j[2] - i[2] * j[1]) + (j[1] - j[2])*x + (i[2] - i[1])*y)) / area_ABC;
                // beta = (0.5f * ((i[2] * j[0] - i[0] * j[2]) + (j[2] - j[0])*x + (i[0] - i[2])*y)) / area_ABC;
                // gamma = (0.5f * ((i[2] * j[0] - i[0] * j[2]) + (j[2] - j[0])*x + (i[0] - i[2])*y)) / area_ABC;
            }

            if (alpha >= 0 && beta >= 0 && gamma >= 0){
                state.image_color[x + state.image_width * y] = make_pixel(255,255,255); 
            }
            
        }
    }

    // for(int y = min_y_Triangle; y < max_y_Triangle; y++){
    //     for (int x = min_x_Triangle; x<max_x_Triangle; x++){
    //         //alpha = area point / abd
    //         float alpha = 
    //     }
    // }
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

