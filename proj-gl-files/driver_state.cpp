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

    state.image_depth = new float[width * height] ;

    for (int i = 0; i < width * height; i++)
    {
	    state.image_depth[i] = 100.0;
    }

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
		//std::cout << "triangles = " <<  (state.num_vertices / 3) << std::endl;
		//std::cout << "num_vertices = " <<  state.num_vertices << std::endl;
		//std::cout << "floats_vertex = " <<  state.floats_per_vertex << std::endl;
                /*for (int i = 0; i < 36 ; i++)
		{
			std::cout << "float [" << i << "] = " << state.vertex_data[i] << std::endl;
		}*/


		for (int k = 0 ; k <  (state.num_vertices / 3); k++)
		{
			data_vertex ver_data;
	        	ver_data.data = new float[MAX_FLOATS_PER_VERTEX];
			const data_geometry *g[3];
			data_geometry geo[3];
			

			for (size_t i = 0; i < 3 ; i++)
			{
				//std::cout << "i = " << i << std::endl;

				for (size_t j = 0; j < state.floats_per_vertex; j++)
				{
				//	std::cout << "index = " << ((state.floats_per_vertex * 3) * k) + (state.floats_per_vertex * i) + j << std::endl;
					ver_data.data[j] = state.vertex_data[((state.floats_per_vertex * 3) * k) + (state.floats_per_vertex * i) + j];
				}
				geo[i].data = new float[MAX_FLOATS_PER_VERTEX];
				state.vertex_shader(ver_data,geo[i],state.uniform_data);

				//geo[i].gl_Position /= geo[i].gl_Position[3];

				g[i] = &geo[i];
			}
			clip_triangle(state,g,);
			//delete geo;
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
    
    float alpha, bravo, charlie;
    data_geometry geo[3];
    const data_geometry *g = new data_geometry[3];
    geo.data = new float[MAX_FLOATS_PER_VERTEX];

    switch (face)
    {
	case 1: // -X
		if (in[0]->gl_Position[0] < (-in[0]->gl_Position[3]))
		{	
			bravo = (((-in[0]->gl_Position[3]) - in[0]->gl_Position[0])/(in[1]->gl_Position[0] + in[1]->gl_Position[3] - in[0]->gl_Position[3] - in[0]->gl_Position[0]));
			charlie = (((-in[0]->gl_Position[3]) - in[0]->gl_Position[0])/(in[2]->gl_Position[0] + in[2]->gl_Position[3] - in[0]->gl_Position[3] - in[0]->gl_Position[0]));
			
			geo[0].data[0] = charlie /in[0]->gl_Position[3];
			geo[0].data[1] = in[0]->gl_Position[1] /in[0]->gl_Position[3];
			geo[0].data[2] = in[0]->gl_Position[2] /in[0]->gl_Position[3];

			geo[1].data[0] = bravo /in[1]->gl_Position[3];
			geo[1].data[1] = in[1]->gl_Position[1] /in[1]->gl_Position[3];
			geo[1].data[2] = in[1]->gl_Position[2] /in[1]->gl_Position[3];

			geo[2].data[0] = in[2]->gl_Position[0] /in[2]->gl_Position[3];
			geo[2].data[1] = in[2]->gl_Position[1] /in[2]->gl_Position[3];
			geo[2].data[2] = in[2]->gl_Position[2] /in[2]->gl_Position[3];
		}
		else if (in[1]->gl_Position[0] < (-in[1]->gl_Position[3]))
		{
			alpha = (((-in[1]->gl_Position[3]) - in[1]->gl_Position[0])/(in[0]->gl_Position[0] + in[0]->gl_Position[3] - in[1]->gl_Position[3] - in[1]->gl_Position[0]));
			charlie = (((-in[1]->gl_Position[3]) - in[1]->gl_Position[0])/(in[2]->gl_Position[0] + in[2]->gl_Position[3] - in[1]->gl_Position[3] - in[1]->gl_Position[0]));
		
			geo[0].data[0] = in[0]->gl_Position[0] /in[0]->gl_Position[3];
			geo[0].data[1] = in[0]->gl_Position[1] /in[0]->gl_Position[3];
			geo[0].data[2] = in[0]->gl_Position[2] /in[0]->gl_Position[3];

			geo[1].data[0] = alpha /in[1]->gl_Position[3];
			geo[1].data[1] = in[1]->gl_Position[1] /in[1]->gl_Position[3];
			geo[1].data[2] = in[1]->gl_Position[2] /in[1]->gl_Position[3];

			geo[2].data[0] = charlie /in[2]->gl_Position[3];
			geo[2].data[1] = in[2]->gl_Position[1] /in[2]->gl_Position[3];
			geo[2].data[2] = in[2]->gl_Position[2] /in[2]->gl_Position[3];
		
		}
		else if (in[2]->gl_Position[0] < (-in[2]->gl_Position[3]))
		{	
			alpha = (((-in[2]->gl_Position[3]) - in[2]->gl_Position[0])/(in[0]->gl_Position[0] + in[0]->gl_Position[3] - in[2]->gl_Position[3] - in[2]->gl_Position[0]));
			bravo = (((-in[2]->gl_Position[3]) - in[2]->gl_Position[0])/(in[1]->gl_Position[0] + in[1]->gl_Position[3] - in[2]->gl_Position[3] - in[2]->gl_Position[0]));

			geo[0].data[0] = alpha / in[0]->gl_Position[3];
			geo[0].data[1] = in[0]->gl_Position[1] / in[0]->gl_Position[3];
			geo[0].data[2] = in[0]->gl_Position[2] / in[0]->gl_Position[3];

			geo[1].data[0] = in[1]->gl_Position[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->gl_Position[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->gl_Position[2] / in[1]->gl_Position[3];

			geo[2].data[0] = bravo / in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}

		break;
	case 0:   // +X
		if (in[0]->data[0] <= in[0]->data[3])
		{	
			bravo = ((in[0]->data[3] - in[0]->data[0])/(in[1]->data[0] - in[1]->data[3] + in[0]->data[3] - in[0]->data[0]));
			charlie = ((in[0]->data[3] - in[0]->data[0])/(in[2]->data[0] - in[2]->data[3] + in[0]->data[3] - in[0]->data[0]));
			
			geo[0].data[0] = charlie /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = bravo /in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		}
		else if (in[1]->data[0] <= in[1]->data[3])
		{
			alpha = ((in[1]->data[3] - in[1]->data[0])/(in[0]->data[0] - in[0]->data[3] + in[1]->data[3] - in[1]->data[0]));
			charlie = ((in[1]->data[3] - in[1]->data[0])/(in[2]->data[0] - in[2]->data[3] + in[1]->data[3] - in[1]->data[0]));
		
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = alpha /in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = charlie /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		
		}
		else if (in[2]->data[0] <= in[2]->data[3])
		{	
			alpha = ((in[2]->data[3] - in[2]->data[0])/(in[0]->data[0] - in[0]->data[3] + in[2]->data[3] - in[2]->data[0]));
			bravo = ((in[2]->data[3] - in[2]->data[0])/(in[1]->data[0] - in[1]->data[3] + in[2]->data[3] - in[2]->data[0]));
	
			geo[0].data[0] = alpha / in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] / in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] / in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] / in[1]->data[3];

			geo[2].data[0] = bravo / in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}
		break;
	case 2: // +Y
 		if (in[0]->data[1] <= in[0]->data[3])
		{	
			bravo = ((in[0]->data[3] - in[0]->data[1])/(in[1]->data[1] - in[1]->data[3] + in[0]->data[3] - in[0]->data[1]));
			charlie = ((in[0]->data[3] - in[0]->data[1])/(in[2]->data[1] - in[2]->data[3] + in[0]->data[3] - in[0]->data[1]));
			
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] =  charlie /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3];
			geo[1].data[1] = bravo /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		}
		else if (in[1]->data[1] <= in[1]->data[3])
		{
			alpha = ((in[1]->data[3] - in[1]->data[1])/(in[0]->data[1] - in[0]->data[3] + in[1]->data[3] - in[1]->data[1]));
			charlie = ((in[1]->data[3] - in[1]->data[1])/(in[2]->data[1] - in[2]->data[3] + in[1]->data[3] - in[1]->data[1]));
		
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = alpha /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] =  charlie /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		
		}
		else if (in[2]->data[1] <= in[2]->data[3])
		{	
			alpha = ((in[2]->data[3] - in[2]->data[1])/(in[0]->data[1] - in[0]->data[3] + in[2]->data[3] - in[2]->data[1]));
			bravo = ((in[2]->data[3] - in[2]->data[1])/(in[1]->data[1] - in[1]->data[3] + in[2]->data[3] - in[2]->data[1]));
	
			geo[0].data[0] = in[0]->data[0] / in[0]->data[3]; 
			geo[0].data[1] = alpha / in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] / in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] / in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] / in[2]->data[3]; 
			geo[2].data[1] = bravo / in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}

		break;
	case 3: // -Y
		if (in[0]->data[1] < (-in[0]->data[3]))
		{	
			bravo = (((-in[0]->data[3]) - in[0]->data[1])/(in[1]->data[1] + in[1]->data[3] - in[0]->data[3] - in[0]->data[1]));
			charlie = (((-in[0]->data[3]) - in[0]->data[1])/(in[2]->data[1] + in[2]->data[3] - in[0]->data[3] - in[0]->data[1]));
			
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3]; 
			geo[0].data[1] = charlie /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = bravo /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		}
		else if (in[1]->data[1] < (-in[1]->data[3]))
		{
			alpha = (((-in[1]->data[3]) - in[1]->data[1])/(in[0]->data[1] + in[0]->data[3] - in[1]->data[3] - in[1]->data[1]));
			charlie = (((-in[1]->data[3]) - in[1]->data[1])/(in[2]->data[1] + in[2]->data[3] - in[1]->data[3] - in[1]->data[1]));
		
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = alpha /in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3]; 
			geo[2].data[1] = charlie /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		
		}
		else if (in[2]->data[1] < (-in[2]->data[3]))
		{	
			alpha = (((-in[2]->data[3]) - in[2]->data[1])/(in[0]->data[1] + in[0]->data[3] - in[2]->data[3] - in[2]->data[1]));
			bravo = (((-in[2]->data[3]) - in[2]->data[1])/(in[1]->data[1] + in[1]->data[3] - in[2]->data[3] - in[2]->data[1]));
	
			geo[0].data[0] = in[0]->data[0] / in[0]->data[3]; 
			geo[0].data[1] = alpha / in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] / in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] / in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] =  bravo / in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}


		break;
	case 4: // +Z
 		if (in[0]->data[2] <= in[0]->data[3])
		{	
			bravo = ((in[0]->data[3] - in[0]->data[2])/(in[1]->data[2] - in[1]->data[3] + in[0]->data[3] - in[0]->data[2]));
			charlie = ((in[0]->data[3] - in[0]->data[2])/(in[2]->data[2] - in[2]->data[3] + in[0]->data[3] - in[0]->data[2]));
			
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3]; 
			geo[0].data[2] =  charlie /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3];
			geo[1].data[2] =  bravo /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		}
		else if (in[1]->data[2] <= in[1]->data[3])
		{
			alpha = ((in[1]->data[3] - in[1]->data[2])/(in[0]->data[2] - in[0]->data[3] + in[1]->data[3] - in[1]->data[2]));
			charlie = ((in[1]->data[3] - in[1]->data[2])/(in[2]->data[2] - in[2]->data[3] + in[1]->data[3] - in[1]->data[2]));
		
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3]; 
			geo[1].data[2] = alpha /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3]; 
			geo[2].data[2] =  charlie /in[2]->data[3];
		
		}
		else if (in[2]->data[2] <= in[2]->data[3])
		{	
			alpha = ((in[2]->data[3] - in[2]->data[2])/(in[0]->data[2] - in[0]->data[3] + in[2]->data[3] - in[2]->data[2]));
			bravo = ((in[2]->data[3] - in[2]->data[2])/(in[1]->data[2] - in[1]->data[3] + in[2]->data[3] - in[2]->data[2]));
	
			geo[0].data[0] = in[0]->data[0] / in[0]->data[3]; 
			geo[0].data[1] = in[0]->data[1] / in[0]->data[3];
			geo[0].data[2] =  alpha / in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] / in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] / in[2]->data[3]; 
			geo[2].data[1] = in[2]->data[1] / in[2]->data[3];
			geo[2].data[2] =  bravo / in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}


		break;
	case 5: // -Z
		if (in[0]->data[2] < (-in[0]->data[3]))
		{	
			bravo = (((-in[0]->data[3]) - in[0]->data[2])/(in[1]->data[2] + in[1]->data[3] - in[0]->data[3] - in[0]->data[2]));
			charlie = (((-in[0]->data[3]) - in[0]->data[2])/(in[2]->data[2] + in[2]->data[3] - in[0]->data[3] - in[0]->data[2]));
			
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3]; 
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] =  charlie /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3];
			geo[1].data[2] =  bravo /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] = in[2]->data[2] /in[2]->data[3];
		}
		else if (in[1]->data[2] < (-in[1]->data[3]))
		{
			alpha = (((-in[1]->data[3]) - in[1]->data[2])/(in[0]->data[2] + in[0]->data[3] - in[1]->data[3] - in[1]->data[2]));
			charlie = (((-in[1]->data[3]) - in[1]->data[2])/(in[2]->data[2] + in[2]->data[3] - in[1]->data[3] - in[1]->data[2]));
		
			geo[0].data[0] = in[0]->data[0] /in[0]->data[3];
			geo[0].data[1] = in[0]->data[1] /in[0]->data[3];
			geo[0].data[2] = in[0]->data[2] /in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] /in[1]->data[3]; 
			geo[1].data[1] = in[1]->data[1] /in[1]->data[3]; 
			geo[1].data[2] = alpha /in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3]; 
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3];
			geo[2].data[2] =  charlie /in[2]->data[3];
		
		}
		else if (in[2]->data[2] < (-in[2]->data[3]))
		{	
			alpha = (((-in[2]->data[3]) - in[2]->data[2])/(in[0]->data[2] + in[0]->data[3] - in[2]->data[3] - in[2]->data[2]));
			bravo = (((-in[2]->data[3]) - in[2]->data[2])/(in[1]->data[2] + in[1]->data[3] - in[2]->data[3] - in[2]->data[2]));
	
			geo[0].data[0] = in[0]->data[0] / in[0]->data[3]; 
			geo[0].data[1] = in[0]->data[1] / in[0]->data[3];
			geo[0].data[2] =  alpha / in[0]->data[3];

			geo[1].data[0] = in[1]->data[0] / in[1]->data[3];
			geo[1].data[1] = in[1]->data[1] / in[1]->data[3];
			geo[1].data[2] = in[1]->data[2] / in[1]->data[3];

			geo[2].data[0] = in[2]->data[0] /in[2]->data[3];
			geo[2].data[1] = in[2]->data[1] /in[2]->data[3]; 
			geo[2].data[2] =  bravo / in[2]->data[3];
			
		}

		for (int i = 0; i < 3; i++ )
		{
			for (int j = 3; j < state.floats_per_vertex; j++)
			{
				geo[i].data[j] = in[i]->data[j];
			}
			g[i] = &geo[i];
		}

		
		rasterize_triangle(state,g);
		return;
		break;
    }

    rasterize_triangle(state,g);
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
			
	
			data_fragment fragment;
			fragment.data = new float[MAX_FLOATS_PER_VERTEX];
			data_output output;
			
			for (int i = 0; i < state.floats_per_vertex; i++)
			{
				switch(state.interp_rules[i])
				{
					case interp_type::flat:
						fragment.data[i] = in[0]->data[i];
						break;
					case interp_type::smooth:
						float Aprime,Bprime,Cprime, alpha2,bravo2,gamma2;
						Aprime = ((alpha * in[0]->data[i])/((alpha * in[0]->data[i])+(bravo * in[1]->data[i])+(gamma * in[2]->data[i]) ));
						Bprime = ((bravo * in[1]->data[i])/((alpha * in[0]->data[i])+(bravo * in[1]->data[i])+(gamma * in[2]->data[i]) ));
						Cprime = ((gamma * in[2]->data[i])/((alpha * in[0]->data[i])+(bravo * in[1]->data[i])+(gamma * in[2]->data[i]) ));
						
						alpha2 = ((Aprime * in[0]->data[i])/((Aprime * in[0]->data[i])+(Bprime * in[1]->data[i])+(Cprime * in[2]->data[i]) )); 
						bravo2 = ((Bprime * in[1]->data[i])/((Aprime * in[0]->data[i])+(Bprime * in[1]->data[i])+(Cprime * in[2]->data[i]) )); 
						gamma2 = ((Cprime * in[2]->data[i])/((Aprime * in[0]->data[i])+(Bprime * in[1]->data[i])+(Cprime * in[2]->data[i]) )); 

						fragment.data[i] = ((alpha2 * in[0]->data[i]) + (bravo2 * in[1]->data[i]) + (gamma2 * in[2]->data[i]));

						break;
					case interp_type::noperspective:
						fragment.data[i] = ((alpha * in[0]->data[i]) + (bravo * in[1]->data[i]) + (gamma * in[2]->data[i]));
						break;	
				}
			}
			state.fragment_shader(fragment,output,state.uniform_data);
			
			float depth = ((in[0]->gl_Position[2] / in[0]->gl_Position[3]) * alpha) +  ((in[1]->gl_Position[2] / in[1]->gl_Position[3]) * bravo) +  ((in[2]->gl_Position[2] / in[2]->gl_Position[3]) * gamma);  

			if (depth < state.image_depth[index])
			{
				state.image_depth[index] = depth; 			
				state.image_color[index] = make_pixel((255 * output.output_color[0]),(255 * output.output_color[1]),(255 * output.output_color[2]));
			}
		}
	}
    }
}

