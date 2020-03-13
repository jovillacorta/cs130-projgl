#include "driver_state.h"
#include <cstring>
#include <float.h>

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

    // our depth for z buffering -- initialize to some large number
    state.image_depth= new float[state.image_height * state.image_width]; 

    state.image_color = new pixel[state.image_height * state.image_width];

    for (unsigned int j = 0; j < state.image_height; j++){
    	for (unsigned int i = 0; i < state.image_width; i++){
    		state.image_color[i + state.image_width * j] = make_pixel(0,0,0);
            state.image_depth[i + state.image_width * j] = FLT_MAX;
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
                clip_triangle(state, triangle, 0);
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

//    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    int A = 0;
    int B = 1;
    int C = 2;
  
    if(face==6)
    {
    //         std::cout << "RASTERIZING triangle with vertices:" << std::endl;
    // std::cout << "A: (" << in[A]->gl_Position[0] << ", " << in[A]->gl_Position[1] << ", " << in[A]->gl_Position[2] << ", " << in[A]->gl_Position[3] << ") " << std::endl; 
    // std::cout << "B: (" << in[B]->gl_Position[0] << ", " << in[B]->gl_Position[1] << ", " << in[B]->gl_Position[2] << ", " << in[B]->gl_Position[3] << ") " << std::endl; 
    // std::cout << "C: (" << in[C]->gl_Position[0] << ", " << in[C]->gl_Position[1] << ", " << in[C]->gl_Position[2] << ", " << in[C]->gl_Position[3] << ") " << std::endl;  
    // std::cout << "A Data: " << *(in[A]->data) << std::endl; 
    // std::cout << "B Data: " << *(in[B]->data) << std::endl; 
    // std::cout << "C Data: " << *(in[C]->data) << std::endl; 
        rasterize_triangle(state, in);
        return;
    }

    // gives us vertice to access depending on the plane
    int clipFace = face % 3;
    // if (clipFace == 0)  std::cout << "Clipping X" << std::endl;
    // if (clipFace == 1)  std::cout << "Clipping Y" << std::endl;
    // if (clipFace == 2)  std::cout << "Clipping Z" << std::endl;
    // do we test against +w or -w? negOrPos tells us which one
    // 0 vs 1
    // 0 is positive, one is negative
    int isNegAxis = face % 2; 
    // if (isNegAxis == 0)  std::cout << "Clipping far" << std::endl;
    // if (isNegAxis == 1)  std::cout << "Clipping near" << std::endl;
    //std::cout << isNegAxis << std::endl;
    int numVertInside = 0;  
    // int X = 0;
    // int Y = 1;
    // int Z = 2;
    int W = 3;


    // std::cout << "Face: " << face << std::endl;
    // std::cout << "clipping triangle with vertices:" << std::endl;
    // std::cout << "A: (" << in[A]->gl_Position[0] << ", " << in[A]->gl_Position[1] << ", " << in[A]->gl_Position[2] << ", " << in[A]->gl_Position[3] << ") " << std::endl; 
    // std::cout << "B: (" << in[B]->gl_Position[0] << ", " << in[B]->gl_Position[1] << ", " << in[B]->gl_Position[2] << ", " << in[B]->gl_Position[3] << ") " << std::endl; 
    // std::cout << "C: (" << in[C]->gl_Position[0] << ", " << in[C]->gl_Position[1] << ", " << in[C]->gl_Position[2] << ", " << in[C]->gl_Position[3] << ") " << std::endl;   


    for (unsigned i = 0; i < 3; i++) {
        // sign will always be either -1 or +1 depending on face
        if (isNegAxis == 0) {
            if ((*in)[i].gl_Position[clipFace] <= (*in)[i].gl_Position[W]){
                numVertInside += 1;
            }
        }
        else {

            if ((*in)[i].gl_Position[clipFace] >= (-1.0f * (*in)[i].gl_Position[W])){
                numVertInside += 1;
            }

        }
    }

    // std::cout << numVertInside << std::endl;
    // if all points are outside, then just completely cull triangle, no need to render
    if (numVertInside == 0) {
        return;
    }

    // if 1 is inside, we clip 2 vertices
    if (numVertInside == 1){
        clip2(state, in, face);
    }
        

    // if 2 are inside, we clip 1 vertices
    if (numVertInside == 2){
        clip1(state, in, face);

    }

    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    else if (numVertInside == 3) {
        clip_triangle(state,in,face+1);
    }


}


// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // note i is i coordinate
    //      j is j coordinate
    //      z is k coordinate for z-buffering
    //      index 0,1,2 is vertice A,B,C respectively
    const int A = 0;
    const int B = 1; 
    const int C = 2;
    float i[3], j[3];

    float z[3]; 

    // calculates i&j coordinates of vertices
    for (int q = 0; q < 3; ++q){

        // adj helps with homogeneous coordinates glPosition 
        // for instance x = x / w
        float adjPosition; 

        // indices from Lab 5:  * note pixels are centered on .5 
        // i = ((x-(-1))/2) * w - .5 
        // j = ((y-(-1))/2) * h - .5
        adjPosition = (in[q]->gl_Position[0]) / (in[q]->gl_Position[3]);
        i[q] = ((adjPosition + 1.0f) / 2.0f * state.image_width) - .5f;
        adjPosition = (in[q]->gl_Position[1]) / (in[q]->gl_Position[3]); 
        j[q] = ((adjPosition + 1.0f) / 2.0f * state.image_height) - .5f;  

        adjPosition = (in[q]->gl_Position[2]) / (in[q]->gl_Position[3]);
        z[q] = ((adjPosition + 1.0f) / 2.0f * state.image_height) - .5f;    
    }



    float area_ABC = (i[A]*(j[B]-j[C]) + i[B]*(j[C]-j[A]) + i[C]*(j[A]-j[B])) * .5f;


    float alpha, beta, gamma;
    float k_alpha[3], k_beta[3], k_gamma[3];

    float zeta;
    // calculate K values for each
    // a = k_0 + k_1*i + k_2*j
    k_alpha[0] = (i[B] * j[C] - i[C] * j[B]) * .5f / area_ABC;
    k_alpha[1] = (j[B] - j[C]) * .5f / area_ABC;
    k_alpha[2] = (i[C] - i[B]) * .5f / area_ABC;

    k_beta[0] = (i[C] * j[A] - i[A] * j[C]) * .5f / area_ABC;
    k_beta[1] = (j[C] - j[A]) * .5f / area_ABC;
    k_beta[2] = (i[A] - i[C]) * .5f / area_ABC;

    k_gamma[0] = (i[A] * j[B] - i[B] * j[A]) * .5f / area_ABC;
    k_gamma[1] = (j[A] - j[B]) * .5f / area_ABC;
    k_gamma[2] = (i[B] - i[A]) * .5f / area_ABC;

    //  instead of iterating through every pixel
    //  we calculate bounding pixels and only iterate through those
    float min_x_Triangle = std::min({i[A],i[B],i[C]});
    float min_y_Triangle = std::min({j[A],j[B],j[C]});
    float max_x_Triangle = std::max({i[A],i[B],i[C]});
    float max_y_Triangle = std::max({j[A],j[B],j[C]});

    // our fragment data
    data_fragment fragment;
    fragment.data = new float[MAX_FLOATS_PER_VERTEX];
    data_output pixel_Color;


    // pseudo code: for all x do 
    //               for all y do
    //                      compute alpha, beta, gamma for (x,y)
    //                          if (0<= alpha,beta,gamma<=1)
    //                              then c = alpha*ca + beta*cb + gamma*cc
    //                                 draw pixel (x,y) with color c
    for (unsigned int y = min_y_Triangle; y < max_y_Triangle; y++){
        for (unsigned int x = min_x_Triangle; x < max_x_Triangle; x++){
            
            for (int v=0; v < 3; v++){
                alpha = k_alpha[0] + k_alpha[1] * x + k_alpha[2] * y;
                beta = k_beta[0] + k_beta[1] * x + k_beta[2] * y;
                gamma = k_gamma[0] + k_gamma[1] * x + k_gamma[2] * y;

            }

            zeta = alpha * z[0] + beta * z[1] + gamma * z[2];
            if (zeta < state.image_depth[x + state.image_width * y] && 
                alpha >= 0 && beta >= 0 && gamma >= 0){

                state.image_depth[x + state.image_width * y] = zeta;

               //state.image_color[x + state.image_width * y] = make_pixel(255,255,255); 
               for (int i = 0; i < state.floats_per_vertex; i++){

               //for smooth type


                switch(state.interp_rules[i]){
                    case interp_type::flat:
                        fragment.data[i] = (*in)[0].data[i];
                        break;
                    case interp_type::smooth:{
                        float alpha1, beta1, gamma1;
                        float k = (alpha / (*in)[0].gl_Position[3])
                            + (beta / (*in)[1].gl_Position[3])
                            + (gamma / (*in)[2].gl_Position[3]);
                        alpha1 = alpha / (k* (*in)[0].gl_Position[3]);
                        beta1 = beta / (k* (*in)[1].gl_Position[3]);    
                        gamma1 = gamma / (k* (*in)[2].gl_Position[3]);

                        fragment.data[i] = alpha1 * (*in)[0].data[i] + beta1 * (*in)[1].data[i] + gamma1 * (*in)[2].data[i];                                        
                        break;
                    }
                    case interp_type::noperspective:
                        fragment.data[i] = alpha * (*in)[0].data[i] + beta * (*in)[1].data[i] + gamma * (*in)[2].data[i];
                        break;
                    default:{}

                }
               }

               state.fragment_shader(fragment, pixel_Color, state.uniform_data);
               state.image_color[x + state.image_width * y] = make_pixel(pixel_Color.output_color[0] * 255,
                                                                        pixel_Color.output_color[1] * 255,
                                                                        pixel_Color.output_color[2] * 255);
            }
            
        }
    }

    delete [] fragment.data;
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

void clip1(driver_state& state, const data_geometry* in[3],int face){
        // gives us vertice to access depending on the plane
    // std::cout << "Clipping two of the "
        // std::cout << "Clipping 1 Vertice" << std::endl;
        int clipFace = face % 3;
        // do we test against +w or -w? negOrPos tells us which one
        // 0 vs 1
        // 0 is positive, one is negative
        int isNegAxis = face % 2; 
            int W = 3;

    int A = 0;
    int B = 1;
    int C = 2;
        
        data_geometry clip1[3]; // out in in

        float alpha1, alpha2;

        if (!isNegAxis) {
            if (((*in)[A].gl_Position[clipFace] > (*in)[A].gl_Position[W])){
                clip1[A] = (*in)[A];
                clip1[B] = (*in)[B];
                clip1[C] = (*in)[C];
                // std::cout << "A is outside" << std::endl;
            }

            else if (((*in)[B].gl_Position[clipFace] > (*in)[B].gl_Position[W])){
                clip1[A] = (*in)[B];
                clip1[B] = (*in)[C];
                clip1[C] = (*in)[A];
                // std::cout << "B is outside" << std::endl;
            }

            else {
                clip1[A] = (*in)[C];
                clip1[B] = (*in)[A];
                clip1[C] = (*in)[B];
                // std::cout << "C is outside" << std::endl;
            }

            alpha1 = (clip1[B].gl_Position[W] - clip1[B].gl_Position[clipFace])
                            / (clip1[A].gl_Position[clipFace] - clip1[A].gl_Position[W] + clip1[B].gl_Position[W] - clip1[B].gl_Position[clipFace]);
            alpha2 = (clip1[C].gl_Position[W] - clip1[C].gl_Position[clipFace])
                            / (clip1[A].gl_Position[clipFace] - clip1[A].gl_Position[W] + clip1[C].gl_Position[W] - clip1[C].gl_Position[clipFace]);
        }

        else {
            if (((*in)[A].gl_Position[clipFace] < (-1 * (*in)[A].gl_Position[W]))){
                clip1[A] = (*in)[A];
                clip1[B] = (*in)[B];
                clip1[C] = (*in)[C];
                // std::cout << "A is outside" << std::endl;
            }

            else if (((*in)[B].gl_Position[clipFace] < (-1 * (*in)[B].gl_Position[W]))){
                clip1[A] = (*in)[B];
                clip1[B] = (*in)[C];
                clip1[C] = (*in)[A];
                // std::cout << "B is outside" << std::endl;
            }

            else {
                clip1[A] = (*in)[C];
                clip1[B] = (*in)[A];
                clip1[C] = (*in)[B];
                // std::cout << "C is outside" << std::endl;
            }
            alpha1 = ((-1 * clip1[B].gl_Position[W]) -  clip1[B].gl_Position[clipFace])
                            / (clip1[A].gl_Position[clipFace] + clip1[A].gl_Position[W] - clip1[B].gl_Position[W] - clip1[B].gl_Position[clipFace]);
            alpha2 = ((-1 * clip1[C].gl_Position[W]) -  clip1[C].gl_Position[clipFace])
                            / (clip1[A].gl_Position[clipFace] + clip1[A].gl_Position[W] - clip1[C].gl_Position[W] - clip1[C].gl_Position[clipFace]);
        }

            vec4 newPoint1 = (alpha1 * clip1[A].gl_Position) + ((1-alpha1)* clip1[B].gl_Position);
            vec4 newPoint2 = (alpha2 * clip1[A].gl_Position) + ((1-alpha2)* clip1[C].gl_Position);

            data_geometry* Triangle1 = new data_geometry[3];
            data_geometry* Triangle2 = new data_geometry[3];

            Triangle1[A].gl_Position = newPoint2;
            Triangle1[A].data = new float[MAX_FLOATS_PER_VERTEX];
            Triangle1[B] = clip1[B];
            Triangle1[C] = clip1[C];  

            Triangle2[A].gl_Position = newPoint1;
            Triangle2[A].data = new float[MAX_FLOATS_PER_VERTEX];
            Triangle2[B] = clip1[B];


            for (int i = 0; i < state.floats_per_vertex; i++){
                switch(state.interp_rules[i]){
                    case interp_type::flat:
                        Triangle1[A].data[i] = clip1[A].data[i];
                        Triangle2[A].data[i] = clip1[A].data[i];
                        break;
                    case interp_type::smooth:{
                            Triangle1[A].data[i] = alpha2 * clip1[A].data[i]  + (1.0f - alpha2) * clip1[C].data[i];  
                            Triangle2[A].data[i] = alpha1 * clip1[A].data[i]  + (1.0f - alpha1) * clip1[B].data[i];                                  
                        break;
                    }
                    case interp_type::noperspective:{
                        float nop_alpha1 = (alpha1 * clip1[A].gl_Position[W]) / (alpha1 * clip1[A].gl_Position[W] + (1.0f - alpha1)* clip1[B].gl_Position[W]);
                        float nop_alpha2 = (alpha2 * clip1[A].gl_Position[W]) / (alpha2 * clip1[A].gl_Position[W] + (1.0f - alpha2)* clip1[C].gl_Position[W]);


                            Triangle2[A].data[i] = nop_alpha1 * clip1[A].data[i]  + (1.0f - nop_alpha1) * clip1[B].data[i];   
                            Triangle1[A].data[i] = nop_alpha2 * clip1[A].data[i] + (1.0f - nop_alpha2) * clip1[C].data[i];
                        
                        break;
                    }
                    default:{}

                }
            }

            Triangle2[C] = Triangle1[A];

            const data_geometry* newTriangle1[3] = {&Triangle1[A], &Triangle1[B], &Triangle1[C]};
            const data_geometry* newTriangle2[3] = {&Triangle2[A], &Triangle2[B], &Triangle2[C]};

            clip_triangle(state, newTriangle1, face +1);
            clip_triangle(state, newTriangle2, face +1);

                delete [] Triangle2[A].data;
                delete [] Triangle1[A].data;
                delete [] Triangle2;
                delete [] Triangle1;
}

void clip2(driver_state& state, const data_geometry* in[3],int face){
    // gives us vertice to access depending on the plane
    int clipFace = face % 3;
    // do we test against +w or -w? negOrPos tells us which one
    // 0 vs 1
    // 0 is positive, one is negative
    int isNegAxis = face % 2; 
        int W = 3;

    int A = 0;
    int B = 1;
    int C = 2;
    data_geometry clip2[3]; // in out out

            float alpha1, alpha2;

            if (!isNegAxis) {
                if ((*in)[A].gl_Position[clipFace] <= (*in)[A].gl_Position[W]){
                    clip2[A] = (*in)[A];
                    clip2[B] = (*in)[B];
                    clip2[C] = (*in)[C];
                }

                else if ((*in)[B].gl_Position[clipFace] <= (*in)[B].gl_Position[W]){
                    clip2[A] = (*in)[B];
                    clip2[B] = (*in)[C];
                    clip2[C] = (*in)[A];
                }

                else {
                    clip2[A] = (*in)[C];
                    clip2[B] = (*in)[A];
                    clip2[C] = (*in)[B];
                }

                alpha1 = (clip2[B].gl_Position[W] - clip2[B].gl_Position[clipFace])
                                / (clip2[A].gl_Position[clipFace] - clip2[A].gl_Position[W] + clip2[B].gl_Position[W] - clip2[B].gl_Position[clipFace]);
                alpha2 = (clip2[C].gl_Position[W] - clip2[C].gl_Position[clipFace])
                                / (clip2[A].gl_Position[clipFace] - clip2[A].gl_Position[W] + clip2[C].gl_Position[W] - clip2[C].gl_Position[clipFace]);
            }

            else {
                if ((*in)[A].gl_Position[clipFace] >= (-1 * (*in)[A].gl_Position[W])){
                    clip2[A] = (*in)[A];
                    clip2[B] = (*in)[B];
                    clip2[C] = (*in)[C];
                }

                else if ((*in)[B].gl_Position[clipFace] >= (-1 * (*in)[B].gl_Position[W])){
                    clip2[A] = (*in)[B];
                    clip2[B] = (*in)[C];
                    clip2[C] = (*in)[A];
                }

                else {
                    clip2[A] = (*in)[C];
                    clip2[B] = (*in)[A];
                    clip2[C] = (*in)[B];
                }
                alpha1 = ((-1.0f * clip2[B].gl_Position[W]) -  clip2[B].gl_Position[clipFace])
                                / (clip2[A].gl_Position[clipFace] + clip2[A].gl_Position[W] - clip2[B].gl_Position[W] - clip2[B].gl_Position[clipFace]);
                alpha2 = ((-1.0f * clip2[C].gl_Position[W]) -  clip2[C].gl_Position[clipFace])
                                / (clip2[A].gl_Position[clipFace] + clip2[A].gl_Position[W] - clip2[C].gl_Position[W] - clip2[C].gl_Position[clipFace]);
            }

                vec4 newPoint1 = (alpha1 * clip2[A].gl_Position) + ((1.0f-alpha1)* clip2[B].gl_Position);
                vec4 newPoint2 = (alpha2 * clip2[A].gl_Position) + ((1.0f-alpha2)* clip2[C].gl_Position);

                data_geometry* Triangle1 = new data_geometry[3];

                Triangle1[A] = clip2[A];
                Triangle1[B].gl_Position = newPoint1;
                Triangle1[B].data = new float[MAX_FLOATS_PER_VERTEX];
                Triangle1[C].gl_Position = newPoint2;
                Triangle1[C].data = new float[MAX_FLOATS_PER_VERTEX];

                for (int i = 0; i < state.floats_per_vertex; i++){
                    switch(state.interp_rules[i]){
                        case interp_type::flat:
                            Triangle1[B].data[i] = clip2[A].data[i];
                            Triangle1[C].data[i] = clip2[A].data[i];
                            break;
                        case interp_type::smooth:{
                                Triangle1[B].data[i] = alpha1 * clip2[A].data[i]  + (1.0f - alpha1) * clip2[B].data[i];  
                                Triangle1[C].data[i] = alpha1 * clip2[A].data[i]  + (1.0f - alpha1) * clip2[C].data[i];                                      
                            break;
                        }
                        case interp_type::noperspective:{
                            float nop_alpha1 = (alpha1 * clip2[A].gl_Position[W]) / (alpha1 * clip2[A].gl_Position[W] + (1.0f - alpha1)* clip2[B].gl_Position[W]);
                            float nop_alpha2 = (alpha2 * clip2[A].gl_Position[W]) / (alpha2 * clip2[A].gl_Position[W] + (1.0f - alpha2)* clip2[C].gl_Position[W]);


                                Triangle1[B].data[i] = nop_alpha1 * clip2[A].data[i]  + (1.0f - nop_alpha1) * clip2[B].data[i];  
                                Triangle1[C].data[i] = nop_alpha2 * clip2[A].data[i] + (1.0f - nop_alpha2) * clip2[C].data[i];
                            break;
                        }
                        default:{}

                    }
                }

                const data_geometry* newTriangle1[3] = {&Triangle1[A], &Triangle1[B], &Triangle1[C]};

                clip_triangle(state, newTriangle1, face +1);

                delete [] Triangle1[B].data;
                delete [] Triangle1[C].data;
                delete [] Triangle1;
}
