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
    state.image_height= height;
    state.image_color = 0;
    state.image_depth = 0;
    state.image_color = new pixel[width*height];
    state.image_depth = new float[width*height];
    for(int i = 0; i < width*height; i++)
    {
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = 2;
    }
}

float Area (float ax, float ay, float bx, float by, float cx, float cy) {
    return 0.5 * ((bx * cy - cx * by) - (ax * cy - cx * ay) + (ax * by - bx * ay));
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
    
    const data_geometry *curr_triangle[3];
    data_geometry geo_arr[3];
    data_vertex vert_dat[3];
    
    switch(type){
        
        
        case (render_type::triangle): { //triangle
        
        int x = 0;
        
            for (unsigned int i = 0; i < state.num_vertices / 3; i++) {
            
                for (unsigned int j = 0; j < 3; j++, x += state.floats_per_vertex) {
                
                    vert_dat[j].data = &state.vertex_data[x];
                    geo_arr[j].data = vert_dat[j].data;
                    state.vertex_shader(vert_dat[j], geo_arr[j], state.uniform_data);
                    curr_triangle[j] = &geo_arr[j];
                    
                }
                
                    clip_triangle(state, curr_triangle, 0);
            }
            
            break;
        
        }
        
        case (render_type::indexed): {
        
            for(int i = 0; i < 3*state.num_triangles; i+=3){
            
                for(int j = 0; j < 3; j++){
                
                    vert_dat[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
                    geo_arr[j].data = vert_dat[j].data;
                    state.vertex_shader(vert_dat[j], geo_arr[j], state.uniform_data);
                    curr_triangle[j] = &geo_arr[j];
                    
                }
                
                clip_triangle(state, curr_triangle, 0);
            }
        
        break;}
        
        
        case (render_type::fan): {
        
            for(int i = 0; i < state.num_vertices; i++){
              
                for(int j = 0; j < 3; j++){
                
                  int index = i + j;
                  if(j == 0) { index = 0; }
                      vert_dat[j].data = &state.vertex_data[index * state.floats_per_vertex];
                  geo_arr[j].data = vert_dat[j].data;
                  state.vertex_shader(vert_dat[j], geo_arr[j], state.uniform_data);
                  curr_triangle[j] = &geo_arr[j];
                 }
            clip_triangle(state, curr_triangle, 0);
            }
                
        break;}
        
        
        case (render_type::strip): { 
        
            for (int i = 0; i < state.num_vertices - 2; ++i) {
            
                for (int j = 0; j < 3; ++j) {
                
                    vert_dat[j].data = &state.vertex_data[(i + j) * state.floats_per_vertex];
                    geo_arr[j].data = vert_dat[j].data;
                    state.vertex_shader(vert_dat[j], geo_arr[j], state.uniform_data);
                    curr_triangle[j] = &geo_arr[j];
                }
                
                
                    clip_triangle(state, curr_triangle, 0);
                
            }
        
        break;}
        
        case (render_type::invalid): { 
        break;}
        
    }
    
    
    
}

//temp function because I'm lazy
void create_new_vertex(driver_state& state, data_geometry* curr_triangle, const data_geometry* v_in, const data_geometry* v_out, int plane, bool pos_sign, float* edited_data) {

	float alpha_smooth = 0, nopersp_alpha = 0;

	if (pos_sign)
		alpha_smooth = (v_out->gl_Position[3] - v_out->gl_Position[plane]) / (v_in->gl_Position[plane] - v_in->gl_Position[3] + v_out->gl_Position[3] - v_out->gl_Position[plane]);
	else
		alpha_smooth = (-v_out->gl_Position[3] - v_out->gl_Position[plane]) / (v_in->gl_Position[plane] + v_in->gl_Position[3] - v_out->gl_Position[3] - v_out->gl_Position[plane]);

	curr_triangle->gl_Position = alpha_smooth * v_in->gl_Position + (1 - alpha_smooth) * v_out->gl_Position;

	nopersp_alpha = alpha_smooth * v_in->gl_Position[3] / (alpha_smooth * v_in->gl_Position[3] + (1 - alpha_smooth) * v_out->gl_Position[3]);

	for (int i = 0; i < state.floats_per_vertex; i++) {
 
		if (state.interp_rules[i] == interp_type::flat) {
   
			edited_data[i] = v_in->data[i];
      
		}
		else if (state.interp_rules[i] == interp_type::smooth) {
   
			edited_data[i] = alpha_smooth * v_in->data[i] + (1 - alpha_smooth) * v_out->data[i];

      
		}
		else if (state.interp_rules[i] == interp_type::noperspective) {
   
			edited_data[i] = nopersp_alpha * v_in->data[i] + (1 - nopersp_alpha) * v_out->data[i];
      
		}
   
   
		}


	curr_triangle->data = edited_data;

}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{

  //declarations for calculations
	bool w_pos, a_in = false, b_in = false, c_in = false;
	int curr_plane = face/2;
 
  if (std::pow(-1, face) > 0){
      w_pos = true;
  }
   else { 
       w_pos = false;
   }   
   

	if (face == 6) {
		
		rasterize_triangle(state, in);
		return;
   
	}
 
   else if (face == 1 || face == 3 || face == 5){ //for negative w_pos 
 
     a_in = (in[0]->gl_Position[curr_plane] >= -in[0]->gl_Position[3]) ? true : false;
		b_in = (in[1]->gl_Position[curr_plane] >= -in[1]->gl_Position[3]) ? true : false;
		c_in = (in[2]->gl_Position[curr_plane] >= -in[2]->gl_Position[3]) ? true : false;
 
   }
   else if (face == 0 || face == 2 || face == 4){ //for positive w_pos
 
     a_in = (in[0]->gl_Position[curr_plane] <= in[0]->gl_Position[3]) ? true : false;
		b_in = (in[1]->gl_Position[curr_plane] <= in[1]->gl_Position[3]) ? true : false;
		c_in = (in[2]->gl_Position[curr_plane] <= in[2]->gl_Position[3]) ? true : false;
 
   }


	//temp data_geometry triangles
	data_geometry* triangle1[3];
	data_geometry* triangle2[3];
 
	for (int i = 0; i < 3; i++) { //initializing
 
		triangle1[i] = new data_geometry();
		triangle2[i] = new data_geometry();
   
	}
 
	float* vertex_data1 = new float[MAX_FLOATS_PER_VERTEX];
	float* vertex_data2 = new float[MAX_FLOATS_PER_VERTEX];
 
     const data_geometry** t1 = const_cast<const data_geometry**>(triangle1);
     const data_geometry** t2 = const_cast<const data_geometry**>(triangle2);
 
 //lazy approach below, could be better

	if (a_in && b_in && c_in) { 
		
		clip_triangle(state, in, face + 1);
   
	}
 else if (!a_in && b_in && !c_in) {

		triangle1[0]->gl_Position = in[1]->gl_Position; 
		triangle1[0]->data = in[1]->data;				
		create_new_vertex(state, triangle1[1], in[1], in[2], curr_plane, w_pos, vertex_data1); 
		create_new_vertex(state, triangle1[2], in[1], in[0], curr_plane, w_pos, vertex_data2); 

		t1 = const_cast<const data_geometry**>(triangle1);
		clip_triangle(state, t1, face + 1);
   
	}
	else if (!a_in && b_in && c_in) {

		triangle1[0]->gl_Position = in[1]->gl_Position;	
		triangle1[0]->data = in[1]->data;				

		triangle1[1]->gl_Position = in[2]->gl_Position; 
		triangle1[1]->data = in[2]->data;	

		create_new_vertex(state, triangle1[2], in[1], in[0], curr_plane, w_pos, vertex_data1); 

		triangle2[0]->gl_Position = in[2]->gl_Position; 
		triangle2[0]->data = in[2]->data;	

		create_new_vertex(state, triangle2[1], in[2], in[0], curr_plane, w_pos, vertex_data2);

		triangle2[2]->gl_Position = triangle1[2]->gl_Position;
		triangle2[2]->data = triangle1[2]->data;			

		t1 = const_cast<const data_geometry**>(triangle1);
		t2 = const_cast<const data_geometry**>(triangle2);
		clip_triangle(state, t1, face + 1);
		clip_triangle(state, t2, face + 1);
	}
 else if (!a_in && !b_in && c_in) {

		triangle1[0]->gl_Position = in[2]->gl_Position; 
		triangle1[0]->data = in[2]->data;			

		create_new_vertex(state, triangle1[1], in[2], in[0], curr_plane, w_pos, vertex_data1); 
		create_new_vertex(state, triangle1[2], in[2], in[1], curr_plane, w_pos, vertex_data2); 

		t1 = const_cast<const data_geometry**>(triangle1);
		clip_triangle(state, t1, face + 1);
   
	}
	else if (a_in && !b_in && c_in) {

		triangle1[0]->gl_Position = in[2]->gl_Position;	
		triangle1[0]->data = in[2]->data;				

		triangle1[1]->gl_Position = in[0]->gl_Position; 
		triangle1[1]->data = in[0]->data;				

		create_new_vertex(state, triangle1[2], in[2], in[1], curr_plane, w_pos, vertex_data1); 

		triangle2[0]->gl_Position = in[0]->gl_Position; 
		triangle2[0]->data = in[0]->data;				

		create_new_vertex(state, triangle2[1], in[0], in[1], curr_plane, w_pos, vertex_data2); 

		triangle2[2]->gl_Position = triangle1[2]->gl_Position;  
		triangle2[2]->data = triangle1[2]->data;				

		t1 = const_cast<const data_geometry**>(triangle1);
		t2 = const_cast<const data_geometry**>(triangle2);
		clip_triangle(state, t1, face + 1);
		clip_triangle(state, t2, face + 1);
   
	}
	else if (a_in && b_in && !c_in) {

		triangle1[0]->gl_Position = in[0]->gl_Position; 
		triangle1[0]->data = in[0]->data;				

		triangle1[1]->gl_Position = in[1]->gl_Position; 
		triangle1[1]->data = in[1]->data;				

		create_new_vertex(state, triangle1[2], in[0], in[2], curr_plane, w_pos, vertex_data1); 

		triangle2[0]->gl_Position = in[1]->gl_Position; 
		triangle2[0]->data = in[1]->data;				

		create_new_vertex(state, triangle2[1], in[1], in[2], curr_plane, w_pos, vertex_data2); 

		triangle2[2]->gl_Position = triangle1[2]->gl_Position;  
		triangle2[2]->data = triangle1[2]->data;				

		t1 = const_cast<const data_geometry**>(triangle1);
		t2 = const_cast<const data_geometry**>(triangle2);
		clip_triangle(state, t1, face + 1);
		clip_triangle(state, t2, face + 1);
   
	}
	else if (a_in && !b_in && !c_in) {

		triangle1[0]->gl_Position = in[0]->gl_Position; 
		triangle1[0]->data = in[0]->data;				
		create_new_vertex(state, triangle1[1], in[0], in[1], curr_plane, w_pos, vertex_data1); 
		create_new_vertex(state, triangle1[2], in[0], in[2], curr_plane, w_pos, vertex_data2); 

		t1 = const_cast<const data_geometry**>(triangle1);
		clip_triangle(state, t1, face + 1);
   
	}


}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    
    
    int x_min, x_max, y_min, y_max; 

    x_min = y_min = 0; // default values of image size

    x_max = state.image_width; 
    y_max = state.image_height;

    float a_x, a_y, b_x, b_y, c_x, c_y;
    float alpha, beta, gamma; 
    vec2 A, B, C; 
    
    a_x = 0.5 * (in[0] -> gl_Position[0] / in[0] -> gl_Position[3] + 1) * state.image_width - 0.5;
    a_y = 0.5 * (in[0] -> gl_Position[1] / in[0] -> gl_Position[3] + 1) * state.image_height - 0.5;
    b_x = 0.5 * (in[1] -> gl_Position[0] / in[1] -> gl_Position[3] + 1) * state.image_width - 0.5;
    b_y = 0.5 * (in[1] -> gl_Position[1] / in[1] -> gl_Position[3] + 1) * state.image_height - 0.5;
    c_x = 0.5 * (in[2] -> gl_Position[0] / in[2] -> gl_Position[3] + 1) * state.image_width - 0.5;
    c_y = 0.5 * (in[2] -> gl_Position[1] / in[2] -> gl_Position[3] + 1) * state.image_height - 0.5;
    
    x_min = std::min(a_x, std::min(b_x, c_x));
    y_min = std::min(a_y, std::min(b_y, c_y));
    x_max = std::max(a_x, std::max(b_x, c_x)) + 1;
    y_max = std::max(a_y, std::max(b_y, c_y)) + 1;
        
    A = {a_x, a_y};
    B = {b_x, b_y};
    C = {c_x, c_y};
    
    float* interpolation_color_data = new float[MAX_FLOATS_PER_VERTEX];

    for (int i = x_min; i < x_max; i++) {
    
        for (int j = y_min; j < y_max; j++) {
            
            vec2 pointVec = {float(i), float(j)};
            
            float abc_area = Area(A[0], A[1], B[0], B[1], C[0], C[1]);

            alpha = Area(pointVec[0], pointVec[1], B[0], B[1], C[0], C[1]) / abc_area;
            beta  = Area(A[0], A[1], pointVec[0], pointVec[1], C[0], C[1]) / abc_area;
            gamma = Area(A[0], A[1], B[0], B[1], pointVec[0], pointVec[1]) / abc_area; 
            
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
            
                float z_point = (alpha * in[0]->gl_Position[2] / in[0]->gl_Position[3]) +(beta * in[1]->gl_Position[2] / in[1]->gl_Position[3]) + (gamma * in[2]->gl_Position[2] / in[2]->gl_Position[3]);
            
                if (z_point < state.image_depth[(state.image_width * j) + i]) {
                
                    state.image_depth[(state.image_width * j) + i] = z_point;
               
                       data_output* return_color = new data_output();
					              data_fragment* color_data = new data_fragment();
            
                        for(int k = 0; k < state.floats_per_vertex; k++){ //interp rule check
                
                          if (state.interp_rules[k] == interp_type::noperspective)
                              interpolation_color_data[k] = alpha*in[0]->data[k] + beta*in[1]->data[k] + gamma*in[2]->data[k];
                          else if (state.interp_rules[k] == interp_type::flat){
                  
                              interpolation_color_data[k] = in[0]->data[k];
                  
                          }
                          else if (state.interp_rules[k] == interp_type::smooth){
                          
                              float alpha_p = 0, beta_p = 0, gamma_p = 0, c = 0;

							                c = (alpha / in[0]->gl_Position[3]) + (beta / in[1]->gl_Position[3]) + (gamma / in[2]->gl_Position[3]);

							                //Calculate the smooth barycentric weights
							                alpha_p = alpha / (in[0]->gl_Position[3] * c);
							                beta_p = beta / (in[1]->gl_Position[3] * c);
							                gamma_p = gamma / (in[2]->gl_Position[3] * c);

							                interpolation_color_data[k] = (alpha_p * in[0]->data[k] + beta_p * in[1]->data[k] + gamma_p * in[2]->data[k]);
                          
                          }
                
                        }
                        color_data->data = interpolation_color_data;

  					      const data_fragment* const_color_data = const_cast<const data_fragment*>(color_data);
					        //Send the interpolated color data to the fragment shader
					        state.fragment_shader(*const_color_data, *return_color, state.uniform_data);

					        //Set the pixel color to the final color
					        state.image_color[(state.image_width * j) + i] =
						        make_pixel(return_color->output_color[0] * 255, return_color->output_color[1] * 255, return_color->output_color[2] * 255);
                        //state.image_color[(j * state.image_width) + i] = make_pixel(255, 255, 255); // if inside the object, color pixel white (for now)
                    
                }
            }
        
          }
        }
}


