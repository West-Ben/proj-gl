#include "driver_state.h"
#include <cstring>
#include <stdlib.h>

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
   
    state.image_color = new pixel[width * height];
    
    
    for (size_t i = 0; i < width * height ; i++)
    {    	
	    state.image_color[i] = make_pixel(0,0,0);
    }

    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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
    //std::cout<< "num_triangles = " << state.num_triangles <<std::endl;

    //std::cout<<"TODO: implement rendering."<<std::endl;
    int count = 0;
    switch (type)
    {
	case render_type::triangle:
		{
//		std::cout << "triangles = " <<  (state.num_vertices / 3) << std::endl;
//		std::cout << "num_vertices = " <<  state.num_vertices << std::endl;
//		std::cout << "floats_vertex = " <<  state.floats_per_vertex << std::endl;
                /*for (int i = 0; i < 36 ; i++)
		{
			std::cout << "float [" << i << "] = " << state.vertex_data[i] << std::endl;
		}*/

		//count++;

		for (int k = 0 ; k <  (state.num_vertices / 3); k++)
		{
			data_vertex ver_data;
	        	ver_data.data = new float[MAX_FLOATS_PER_VERTEX];
			const data_geometry *g[3];
			data_geometry geo[3];
			for (size_t i = 0; i < 3 ; i++)
			{
				for (size_t j = 0; j < state.floats_per_vertex; j++)
				{
					//std::cout << "index = " << ((state.floats_per_vertex * 3) * k) + (state.floats_per_vertex * i) + j << std::endl;
					ver_data.data[j] = state.vertex_data[((state.floats_per_vertex * 3) * k) + (state.floats_per_vertex * i) + j];
				}
				state.vertex_shader(ver_data,geo[i],state.uniform_data);
				geo[i].gl_Position /= geo[i].gl_Position[3];

				g[i] = &geo[i];
			}
			rasterize_triangle(state,g);
		}
		break;	
    		}
	case render_type::indexed:
		{
		break;
    		}
	case render_type::fan:
		{
		//std::cout << "fan count = " << count << std::endl;
		break;
    		}	
	case render_type::strip:
		{
		//std::cout << "strip count = " << count << std::endl;
		break;
    			
		}
    }

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
//    std::cout<<"TODO: implement rasterization"<<std::endl;

   float ai = (((state.image_width/2) * in[0]->gl_Position[0]) + ((state.image_width/2) - 0.5));
   float aj = (((state.image_height/2) * in[0]->gl_Position[1]) + ((state.image_height/2) - 0.5));

   float bi =  (((state.image_width/2) * in[1]->gl_Position[0]) + ((state.image_width/2) - 0.5));
   float bj = (((state.image_height/2) * in[1]->gl_Position[1]) + ((state.image_height/2) - 0.5));

   float ci =  (((state.image_width/2) * in[2]->gl_Position[0]) + ((state.image_width/2) - 0.5));
   float cj = (((state.image_height/2) * in[2]->gl_Position[1]) + ((state.image_height/2) - 0.5));

   float area = 0.5 * (((bi * cj) - (ci * bj)) - 
			((ai * cj) - (ci * aj)) +  
			((ai * bj) - (bi * aj))); 

    int xmin = std::min(std::min(ai,bi),ci);
    int xmax = std::max(std::max(ai,bi),ci);
    int ymin = std::min(std::min(aj,bj),cj);
    int ymax = std::max(std::max(aj,bj),cj);
    for (int i = xmin; i < xmax; i++)
    {

    	for (int j = ymin; j < ymax ; j++)
	{
    		float alpha = 0.5 * (((bi * cj) - (ci * bj)) - 
			((i * cj) - (ci * j)) +  
			((i * bj) - (bi * j))); 
    		float bravo = 0.5 * (((i * cj) - (ci * j)) - 
			((ai * cj) - (ci * aj)) +  
			((ai * j) - (i * aj))); 
    		float gamma = 0.5 * (((bi * j) - (i * bj)) - 
			((ai * j) - (i * aj)) +  
			((ai * bj) - (bi * aj))); 


		alpha /= area;
		bravo /= area;
		gamma /= area;
		
		if (alpha >= 0 && bravo >= 0 && gamma >= 0 && gamma >= 0 && alpha + bravo + gamma <= 1.01)
		{

			int index = (state.image_width  * j) + i;
			state.image_color[index] = make_pixel(255,255,255);
		}
	}
    }


}

