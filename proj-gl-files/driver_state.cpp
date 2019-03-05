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
			
			//std::cout<<"geo[0]->gl_Position[3] = " << geo[0].gl_Position[3]<<std::endl;
			//std::cout<<"geo[1]->gl_Position[3] = " << geo[1].gl_Position[3]<<std::endl;
			//std::cout<<"geo[2]->gl_Position[3] = " << geo[2].gl_Position[3]<<std::endl;
			clip_triangle(state,g,0);
			//delete geo;
		}
		break;	
    	}
	case render_type::indexed:
		{
			for (int k = 0 ; k <  state.num_triangles; k++)
			{
				data_vertex ver_data;
				ver_data.data = new float[MAX_FLOATS_PER_VERTEX];
				
				const data_geometry *g[3];
				data_geometry geo[3];
				

				for (size_t i = 0; i < 3 ; i++)
				{
					for(int j = 0; j < state.floats_per_vertex; j++)
					{
						ver_data.data[j] = state.vertex_data[state.index_data[(k * 3) + i] + j];
					}

					
					geo[i].data = new float[MAX_FLOATS_PER_VERTEX];
					state.vertex_shader(ver_data,geo[i],state.uniform_data);
					
					g[i] = &geo[i];
				}
				
				clip_triangle(state,g,0);
			}
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
	data_geometry geo[3];
    const data_geometry *g[3];
	
    if(face==6)
    {
		for (int i = 0; i < 3; i++)
		{
			geo[i].data = in[i]->data;
			geo[i].gl_Position = (in[i]->gl_Position / in[i]->gl_Position[3]);
			
			g[i] = &geo[i];
		}
		
        rasterize_triangle(state, g);
        return;
    }

	//std::cout<<"Sin[0]->gl_Position[3] = " << in[0]->gl_Position[3]<<std::endl;
	//std::cout<<"Sin[1]->gl_Position[3] = " << in[1]->gl_Position[3]<<std::endl;
	//std::cout<<"Sin[2]->gl_Position[3] = " << in[2]->gl_Position[3]<<std::endl;
    float alpha, bravo, charlie;
	vec4 Pa, Pb, Pc;

    int axis;
	int sign;
	
	for (int i = 0; i < 3; i++ )
	{
		geo[i].data = in[i]->data;
	}


	switch (face)
	{
		case 0:
			axis = 0;
			sign = 1;
		break;
		case 1:
			axis = 0;
			sign = -1;
		break;
		case 2:
			axis = 1;
			sign = 1;
		break;
		case 3:
			axis = 1;
			sign = -1;
		break;
		case 4:
			axis = 2;
			sign = 1;
		break;
		case 5:
			axis = 2;
			sign = -1;
		break;
	}

	int abc[3];

	if (sign == 1)
	{
		for (int vert = 0; vert < 3; vert++)
		{
			//std::cout << " in["<< vert <<"]->gl_Position["<< axis <<"] = " << in[vert]->gl_Position[axis] << std::endl;
			if (in[vert]->gl_Position[axis] > in[vert]->gl_Position[3])
			{
				abc[vert] = 0;
			}
			else
			{
				abc[vert] = 1;
			}
		}
	}
	else
	{
		for (int vert = 0; vert < 3; vert++)
		{
			//std::cout << " - in["<< vert <<"]->gl_Position["<< axis <<"] = " << in[vert]->gl_Position[axis] << std::endl;
			if (in[vert]->gl_Position[axis] < -in[vert]->gl_Position[3])
			{
				abc[vert] = 0;
			}
			else
			{
				abc[vert] = 1;
			}
		}
	}

	if( abc[0] == 0 && abc[1] == 0 && abc[2] == 1)
	{
		alpha = (((sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis])/(in[2]->gl_Position[axis] - (sign * in[2]->gl_Position[3]) + (sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis]));
		bravo = (((sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis])/(in[2]->gl_Position[axis] - (sign * in[2]->gl_Position[3]) + (sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis]));
		
		Pa = (alpha * in[2]->gl_Position) + ((1 - alpha) * in[0]->gl_Position);
		Pb = (bravo * in[2]->gl_Position) + ((1 - bravo) * in[1]->gl_Position);
		
		//std::cout << "pa = " << Pa << std::endl;
		//std::cout << "pb = " << Pb << std::endl;
		
		geo[0].gl_Position = Pa;
		geo[1].gl_Position = Pb;
		geo[2].gl_Position = in[2]->gl_Position;

		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "-X" << std::endl;
		clip_triangle(state,g,face+1);
		return;
	}
	else if(abc[0] == 0 && abc[1] == 1 && abc[2] == 0)
	{
		alpha = (((sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis])/(in[1]->gl_Position[axis] - (sign * in[1]->gl_Position[3]) + (sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis]));
		charlie = (((sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis])/(in[1]->gl_Position[axis] - (sign * in[1]->gl_Position[3]) + (sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis]));
		
		Pa = (alpha * in[1]->gl_Position) + ((1 - alpha) * in[0]->gl_Position);
		Pc = (charlie * in[1]->gl_Position) + ((1 - charlie) * in[2]->gl_Position);
		

		//std::cout << "pa = " << Pa << std::endl;
		//std::cout << "pc = " << Pc << std::endl;
		
		geo[0].gl_Position = Pa;
		geo[1].gl_Position = in[1]->gl_Position;
		geo[2].gl_Position = Pc;

		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "b in" << std::endl;
		clip_triangle(state,g,face+1);
		return;
	}
	else if (abc[0] == 0 && abc[1] == 1 && abc[2] == 1)
	{
		bravo = (((sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis])/(in[1]->gl_Position[axis] - (sign * in[1]->gl_Position[3]) + (sign * in[0]->gl_Position[3]) - in[0]->gl_Position[axis]));
		charlie = (((sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis])/(in[0]->gl_Position[axis] - (sign * in[0]->gl_Position[3]) + (sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis]));
		
		Pb = (bravo * in[1]->gl_Position) + ((1 - bravo) * in[0]->gl_Position);
		Pc = (charlie * in[0]->gl_Position) + ((1 - charlie) * in[2]->gl_Position);
		
		geo[0].gl_Position = Pc;
		geo[1].gl_Position = in[1]->gl_Position;
		geo[2].gl_Position = Pb;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "b and c in" << std::endl;
		clip_triangle(state,g,face+1);
		
		geo[0].gl_Position = Pc;
		geo[1].gl_Position = in[1]->gl_Position;
		geo[2].gl_Position = in[2]->gl_Position;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "2 b and c in" << std::endl;
		clip_triangle(state,g,face+1);
		
		return;
	}
	else if (abc[0] == 1 && abc[1] == 0 && abc[2] == 0)
	{
		bravo = (((sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis])/(in[0]->gl_Position[axis] - (sign * in[0]->gl_Position[3]) + (sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis]));
		charlie = (((sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis])/(in[0]->gl_Position[axis] - (sign * in[0]->gl_Position[3]) + (sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis]));
		
		Pb = (bravo * in[0]->gl_Position) + ((1 - bravo) * in[1]->gl_Position);
		Pc = (charlie * in[0]->gl_Position) + ((1 - charlie) * in[2]->gl_Position);
		
		geo[0].gl_Position = in[0]->gl_Position;
		geo[1].gl_Position = Pb;
		geo[2].gl_Position = Pc;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "a in" << std::endl;
		clip_triangle(state,g,face+1);
		return;
	}
	else if (abc[0] == 1 && abc[1] == 0 && abc[2] == 1)
	{
		alpha = (((sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis])/(in[0]->gl_Position[axis] - (sign * in[0]->gl_Position[3]) + (sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis]));
		charlie = (((sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis])/(in[2]->gl_Position[axis] - (sign * in[2]->gl_Position[3]) + (sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis]));
		
		Pa = (alpha * in[0]->gl_Position) + ((1 - alpha) * in[1]->gl_Position);
		Pc = (charlie * in[2]->gl_Position) + ((1 - charlie) * in[1]->gl_Position);
		
		geo[0].gl_Position = in[0]->gl_Position;
		geo[1].gl_Position = Pa;
		geo[2].gl_Position = Pc;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "a and c in" << std::endl;
		clip_triangle(state,g,face+1);
		
		geo[0].gl_Position = in[0]->gl_Position;
		geo[1].gl_Position = Pc;
		geo[2].gl_Position = in[2]->gl_Position;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "2 a and c in" << std::endl;
		clip_triangle(state,g,face+1);
		return;
	}
	else if (abc[0] == 1 && abc[1] == 1 && abc[2] == 0)
	{
		alpha = (((sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis])/(in[0]->gl_Position[axis] - (sign * in[0]->gl_Position[3]) + (sign * in[2]->gl_Position[3]) - in[2]->gl_Position[axis]));
		bravo = (((sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis])/(in[2]->gl_Position[axis] - (sign * in[2]->gl_Position[3]) + (sign * in[1]->gl_Position[3]) - in[1]->gl_Position[axis]));
		
		Pa = (alpha * in[0]->gl_Position) + ((1 - alpha) * in[2]->gl_Position);
		Pb = (bravo * in[2]->gl_Position) + ((1 - bravo) * in[1]->gl_Position);
		
		geo[0].gl_Position = in[0]->gl_Position;
		geo[1].gl_Position = Pa;
		geo[2].gl_Position = Pb;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "a and b in" << std::endl;
		clip_triangle(state,g,face+1);
		
		geo[0].gl_Position = in[0]->gl_Position;
		geo[1].gl_Position = in[1]->gl_Position;
		geo[2].gl_Position = Pb;
		
		for (int i = 0; i < 3; i++ )
		{
			g[i] = &geo[i];
			
		}
		
		//std::cout << "2 a and b in" << std::endl;
		clip_triangle(state,g,face+1);
		return;
	}
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
//    std::cout<<"TODO: implement rasterization"<<std::endl;
/*for(int i = 0; i < 3;i++)
{
	for(int j=0; j < 3;j++)
	{
		std::cout<<"in["<< i <<"]->gl_Position["<< j <<"] = "<< in[i]->gl_Position[j] <<std::endl;
	}
}*/
	
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
	
	//std::cout<<"xmin = "<< xmin <<std::endl;
	//std::cout<<"xmax = "<< xmax <<std::endl;
	//std::cout<<"ymin = "<< ymin <<std::endl;
	//std::cout<<"ymax = "<< ymax <<std::endl;
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

			//std::cout<<"alpha + bravo + gamma = "<< alpha + bravo + gamma <<std::endl;

			alpha /= area;
			bravo /= area;
			gamma /= area;
			//std::cout<<"alpha + bravo + gamma = "<< alpha + bravo + gamma <<std::endl;
			
			if (alpha >= 0 && bravo >= 0 && gamma >= 0 && gamma >= 0 && alpha + bravo + gamma <= 1.0001 )
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

