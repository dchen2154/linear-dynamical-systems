/**
 * File: linear_dynamical_systems.pde
 * Project: Visualizing Linear Dynamical Systems
 *
 * Given a linear dynamical system of form
 *    x_{n+1} = A x_{n}
 * by having the user specify the matrix A and the initial state vector x_{0},
 * this app draws a diagram to depict the states x_{n} as the number of
 * iterations n increase. The matrix and initial state can be
 * changed in the configuration mode. There, you also have the option to change
 * some visual settings. After the user sets their desired configuration and
 * continues, the app depicts the diagram drawn as a graph with vertices
 * corresponding to each column of the matrix. As the number of iterations
 * increase, an animation will play to show the proportional change in states
 * and how a component affects the component of the next state. The app also
 * keeps a local history to return to previous iterations (but not systems).
 * When rewinding, no animations will play.
 *
 * Author: David Chen <dchen2@andrew.cmu.edu>
 */


/********** Helper Functions **********/

/**
 * Returns the product of two matrices or throws an error if matrices have
 * incompatible size.
 */
float[][] matrixMultiply(float[][] A, float[][] B) {
  if (A.length != B[0].length)
    throw new RuntimeException("Attempt to multiply matrices of invalid sizes");
  float[][] M = new float[B.length][A[0].length];
  for (int i = 0; i < B.length; i++) {
    for (int j = 0; j < A[0].length; j++) {
      float entry = 0;
      for (int k = 0; k < A.length; k++)
        entry += A[k][j] * B[i][k];
      M[i][j] = entry;
    }
  }
  return M;
}

/**
 * Given a vector and an adjusted total, calculate the ratio of the adjusted
 * total to the actual total. Essentially, this is the scale factor for which
 * the sum of the elements of the scaled vector equals the adjusted total.
 */
float weightRatio(float[][] v, float adjTotal) {
  if (v.length != 1)
    throw new RuntimeException("Invalid vector");
  float actTotal = 0;
  for (int i = 0; i < v[0].length; i++)
    actTotal += v[0][i];
  return adjTotal / actTotal;
}

/**
 * Checks if x is between b1 and b2.
 */
boolean inbtw(float b1, float b2, float x) {
  return b1 <= x && x <= b2 || b2 <= x && x <= b1;
}

/**
 * Returns the string with c inserted at index i of string s
 */
String string_append(String s, int i, char c) {
  if (i < 0 || i > s.length())
    throw new ArrayIndexOutOfBoundsException();
  return s.substring(0, i) + c + s.substring(i);
}

/**
 * Returns the string with index i removed from string s
 */
String string_remove(String s, int i) {
  if (i < 0 || i >= s.length())
    throw new ArrayIndexOutOfBoundsException();
  return s.substring(0, i) + s.substring(i + 1);
}


/********** General Variables **********/

/**
 * Enum type for current mode of the app
 *    Menu
 *    Help
 *    Configuration
 *    Main (where the animation occurs)
 */
app_mode mode;

/** Text size for every mode except main */
int text_size;

/** Spacing values for text */
float text_margin, text_indent;


/********** Help Mode Variables **********/

/** Help message header and body */
String help_hdr, help_msg;


/********** Configuration Mode Variables **********/

/**
 * Enum type for current field of the configuration option
 *
 * In the configuration mode, the user can change the follow properties.
 *    Dimension
 *    Matrix
 *    Vector
 *    Time Step Size
 *    Total Nodes
 *    Node Spread Ratio
 *    Continue
 * The Continue field acts as a jump to the animation; it is not an actual
 * editable field
 */
config_field config_editor_field;

/** Current column of the editor */
int config_editor_entry;

/** Scroll values */
float config_scroll_x, config_scroll_y;

/** String representation of the dimension count of the system */
String config_dim_str;

/** Location of the dimensions value */
float config_dim_x, config_dim_y;

/** String representation of the matrix of the system */
String[][] config_mat_str;

/** Location of the matrix values */
float config_mat_x, config_mat_y;

/**
 * List of coordinates for each block of the matrix
 *
 * config_mat_coords[i][j][0] = x-coordinate of the (i, j)th entry
 * config_mat_coords[i][j][1] = y-coordinate of the (i, j)th entry
 */
float[][][] config_mat_coords;

/** Current index of the (string) matrix */
int[] config_mat_cind;

/**
 * Width of the corresponding vector of a matrix entry
 *
 * config_mat_w[i] = width of the i-th column of config_mat_str
 */
float[] config_mat_w;

/** String representation of the initial vector of the system */
String[] config_vec_str;

/** Location of the initial state vector values */
float config_vec_x, config_vec_y;

/**
 * List of coordinates for each block of the vector
 *
 * config_vec_x         = x-coordinate of every entry
 * config_vec_coords[i] = y-coordinate of the i-th entry
 */
float[] config_vec_coords;

/** Current index of the (string) vector */
int config_vec_cind;

/** String representation of the time step size */
String config_dt_str;

/** Location of the time step size value */
float config_dt_x, config_dt_y;

/** String representation of the total nodes */
String config_ntotal_str;

/** Location of the total nodes value */
float config_ntotal_x, config_ntotal_y;

/** String representation of the node spread ratio */
String config_nspread_str;

/** Location of the node spread ratio value */
float config_nspread_x, config_nspread_y;

/** Location of the continue option */
float config_cont_x, config_cont_y;


/********** User-Inputted Variables **********/

/**
 * Number of dimensions of the system
 *
 * Invariant: state_dim == state_mat.length == state.mat[0].length
 */
int state_dim;

/** Matrix and initial vector of the system */
float[][] state_mat, state_vec;

/** Time step size */
float state_dt;

/**
 * Total number of nodes
 *
 * Initial number of nodes in the system. If the sum of the columns of the
 * matrix equal 1 then this is also the total number of nodes. Note that the
 * animation rounds up, so the user might see more.
 */
int state_ntotal;

/** 
 * Ratio of total distance for which nodes can occupy
 *
 * Intuitively, this is a node spread ratio
 */
float state_nspread;


/********** Main Mode Variables **********/

/**
 * Animation time
 *
 * At t = 1, the first nodes sent will reach their destination.
 * At t = 1 + state_nspread, the last nodes sent will reach their destination.
 */
float t;

/** Iteration count */
int iters;

/** Center of the canvas */
float cx, cy;

/** Distance from the center for which vertices are drawn */
float cr;

/** Distance from the center for which the values at each vertices are drawn */
float tr;

/** Text size for main mode */
int vert_textsize;

/** Diameter of each vertex and diameter of each node */
float vert_d, node_d;

/**
 * For each vertex, get the position of itself and its value
 *
 * vert_pos[i][0]  = x-coordinate of the i-th vertex
 * vert_pos[i][1]  = y-coordinate of the i-th vertex
 * vert_tpos[i][0] = x-coordinate of the i-th vertex's value
 * vert_tpos[i][1] = y-coordinate of the i-th vertex's value
 */
float[][] vert_pos, vert_tpos;

/**
 * For each directed edge, appoint a color
 *
 * vert_colors[i][j] = color from vertex i to vertex j
 */
float[][] vert_colors;

/**
 * For each directed edge, note the intial and final intial
 *
 * node_init[i][0] = x-coordinate of the start of the i-th edge
 * node_init[i][1] = y-coordinate of the start of the i-th edge
 * node_finl[i][0] = x-coordinate of the end   of the i-th edge
 * node_finl[i][1] = y-coordinate of the end   of the i-th edge
 */
float[][] node_init, node_finl;

/**
 * For each corresponding directed edge, record the number of nodes traveling
 *
 * node_counts[i] = number of nodes along edge i
 */
float[] node_counts;

/** Length of node_init, node_finl, node_counts */
int node_len;

/** Boolean to lock all events whilst animating */
boolean isAnimating;


/********** History Variables **********/

/** History of node_init, node_finl, and the states respectively */
ArrayList<float[][]> hist_ninit, hist_nfinl, hist_states;

/** History of node_counts */
ArrayList<float[]> hist_ncnts;

/**
 * History of node_len
 *
 * hist_ncnts[i] = length(hist_ninit[i])
 *               = length(hist_nfinl[i])
 *               = length(hist_ncnts[i])
 */
IntList hist_nlen;


/********** Setup Functions **********/

/**
 * Sets up the help mode
 *
 * Called only upon launching to grab the help strings.
 */
void setup_help() {
  help_hdr = "Welcome to Visualizing Linear Dynamic Systems!";
  help_msg = loadStrings("help_msg.txt")[0];
}

/**
 * Sets up the configurations mode
 *
 * Called only upon launching to set the config variables.
 */
void setup_config() {
  config_editor_field = config_field.DIM;
  config_editor_entry = 0;
  config_scroll_x     = 0;
  config_scroll_y     = 0;

  config_dim_str      = "";
  config_dim_x        = text_indent + textWidth("Dimension:") + text_size;
  config_dim_y        = text_margin + 2 * text_size;

  config_mat_str      = new String[0][0];
  config_mat_x        = text_indent + textWidth("Matrix:") + text_size;
  config_mat_y        = text_margin + 4 * text_size;
  config_mat_coords   = new float[0][0][0];
  config_mat_cind     = new int[2];
  config_mat_cind[0]  = 0;
  config_mat_cind[1]  = 0;
  config_mat_w        = new float[0];

  config_vec_str      = new String[0];
  config_vec_x        = text_indent + textWidth("Initial State:") + text_size;
  config_vec_y        = text_margin + 6 * text_size;
  config_vec_coords   = new float[0];
  config_vec_cind     = 0;

  config_dt_str       = "";
  config_dt_x         = text_indent + textWidth("Time Step Size:") + text_size;
  config_dt_y         = text_margin + 8 * text_size;

  config_ntotal_str   = "";
  config_ntotal_x     = text_indent + textWidth("Total Nodes:") + text_size;
  config_ntotal_y     = text_margin + 10 * text_size;

  config_nspread_str  = "";
  config_nspread_x    = text_indent + textWidth("Node Spread Ratio:") + text_size;
  config_nspread_y    = text_margin + 12 * text_size;

  config_cont_x       = text_margin + textWidth("Press here to continue.");
  config_cont_y       = text_margin + 14 * text_size;
}

/**
 * Sets up the main mode
 *
 * Called upon each entrance into main mode.
 */
void setup_main() {
  int n     = state_dim * (state_dim - 1);
  t         = 0;
  iters     = 0;
  vert_pos  = new float[state_dim][2];
  vert_tpos = new float[state_dim][2];
  if (state_dim == 1) {
    vert_d          = 1.5 * text_size;
    vert_pos[0][0]  = cx;
    vert_pos[0][1]  = cy;
    vert_tpos[0][0] = cx;
    vert_tpos[0][1] = cy + vert_d / 2 + vert_textsize;
  } else {
    vert_d = min(1.5 * text_size, 0.5 * 2 * cr * sin(PI / state_dim));
    for (int i = 0; i < state_dim; i++) {
      float ang = i * TAU / state_dim;
      vert_pos[i][0]  = cx + cr * sin(ang);
      vert_pos[i][1]  = cy - cr * cos(ang);
      vert_tpos[i][0] = cx + tr * sin(ang);
      vert_tpos[i][1] = cy - tr * cos(ang);
    }
  }
  node_d        = vert_d / 2;
  vert_textsize = int(2.0 / 3.0 * vert_d);
  vert_colors   = new float[n][3];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 3; j++)
      vert_colors[i][j] = random(0, 256);
  }
  isAnimating = false;
  hist_ninit  = new ArrayList<float[][]>();
  hist_nfinl  = new ArrayList<float[][]>();
  hist_states = new ArrayList<float[][]>();
  hist_ncnts  = new ArrayList<float[]>();
  hist_nlen   = new IntList();
  hist_states.add(state_vec);
}

/** Setup */
void setup() {
  size(800, 800);
  background(0);
  frameRate(60);
  int d = min(width, height);

  mode = app_mode.MENU;

  text_size = d / 25;
  textSize(text_size);
  text_margin = text_size / 2;
  text_indent = text_margin + text_size;

  setup_config();

  setup_help();

  cx = width / 2;
  cy = height / 2;
  cr = 0.3 * min(width, height);
  tr = 0.4 * min(width, height);
  setup_main();
}


/********** Event Functions for Configuration Mode **********/

/**
 * Updates state_dim according to the user input and recreates the other config
 * variables accordingly
 *
 * Called after the user leaves the Dimesion field
 */
void update_dim() {
  if (config_dim_str.equals("") || config_dim_str.equals("0"))
    config_dim_str = "1";
  if (state_dim == int(config_dim_str))
    return;
  state_dim         = int(config_dim_str);
  config_dim_str    = str(state_dim);
  config_mat_str    = new String[state_dim][state_dim];
  config_vec_str    = new String[state_dim];
  config_mat_coords = new float[state_dim][state_dim][2];
  config_vec_coords = new float[state_dim];
  config_mat_w      = new float[state_dim];
  config_vec_y      = text_margin + (4 + 2 * state_dim) * text_size;
  float w = 0;
  for (int i = 0; i < state_dim; i++) {
    config_mat_w[i] = textWidth("_");
    for (int j = 0; j < state_dim; j++) {
      config_mat_str[i][j]       = "";
      config_mat_coords[i][j][0] = config_mat_x + w + i * text_size;
      config_mat_coords[i][j][1] = config_mat_y + 2 * j * text_size;
    }
    config_vec_str[i]    = "";
    config_vec_coords[i] = config_vec_y + 2 * i * text_size;
    w += config_mat_w[i];
  }
  config_dt_y      = config_vec_y + (2 * state_dim) * text_size;
  config_ntotal_y  = config_dt_y + 2 * text_size;
  config_nspread_y = config_ntotal_y + 2 * text_size;
  config_cont_y    = config_nspread_y + 2 * text_size;
}

/**
 * Editor for the Dimension field
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_dim() {
  if (key == CODED) {
    if (keyCode == UP) {
      config_editor_field = config_field.CONT;
      config_editor_entry = 0;
      update_dim();
    } else if (keyCode == DOWN) {
      config_editor_field = config_field.MAT;
      config_editor_entry = 0;
      config_mat_cind[0]  = 0;
      config_mat_cind[1]  = 0;
      update_dim();
    } else if (keyCode == LEFT) {
        if (config_editor_entry == 0) {
          config_editor_field = config_field.CONT;
          config_editor_entry = 0;
          update_dim();
        } else
          config_editor_entry--;
    } else if (keyCode == RIGHT) {
        if (config_editor_entry == config_dim_str.length()) {
          config_editor_field = config_field.MAT;
          config_mat_cind[0]  = 0;
          config_mat_cind[1]  = 0;
          config_editor_entry = 0;
          update_dim();
        } else
          config_editor_entry++;
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      config_editor_field = config_field.MAT;
      update_dim();
      config_mat_cind[0]  = 0;
      config_mat_cind[1]  = 0;
      config_editor_entry = config_mat_str[0][0].length();
    } else if (key == 8 && config_editor_entry > 0) {
      // Backspace
      config_dim_str = string_remove(config_dim_str, config_editor_entry - 1);
      config_editor_entry--;
    } else if (key == 127 && config_editor_entry < config_dim_str.length()) {
      // Delete
      config_dim_str = string_remove(config_dim_str, config_editor_entry);
    } else if ('0' <= key && key <= '9') {
      config_dim_str = string_append(config_dim_str, config_editor_entry, key);
      config_editor_entry++;
    }
  }
}

/**
 * Updates state_mat based on user input
 *
 * Called after the user leaves the Matrix field.
 */
void update_mat() {
  state_mat = new float[state_dim][state_dim];
  for (int i = 0; i < state_dim; i++) {
    for (int j = 0; j < state_dim; j++) {
      if (config_mat_str[i][j].equals("") || config_mat_str[i][j].equals("."))
        config_mat_str[i][j] = "0";
      state_mat[j][i] = float(config_mat_str[i][j]);
      if (state_mat[j][i] == int(state_mat[j][i]))
        config_mat_str[i][j] = str(int(state_mat[j][i]));
      else
        config_mat_str[i][j] = str(state_mat[j][i]);
    }
  }
  update_strmat();
}

/**
 * Updates config_mat_str based on user input
 *
 * Called everytime the matrix is altered in the config mode.
 */
void update_strmat() {
  for (int i = 0; i < state_dim; i++) {
    config_mat_w[i] = textWidth("_");
    for (int j = 0; j < state_dim; j++)
      config_mat_w[i] = max(config_mat_w[i], textWidth(config_mat_str[i][j]));
  }
  float w = 0;
  for (int i = 0; i < state_dim; i++) {
    for (int j = 0; j < state_dim; j++) {
      config_mat_coords[i][j][0] = config_mat_x + w + i * text_size;
      config_mat_coords[i][j][1] = config_mat_y + 2 * j * text_size;
    }
    w += config_mat_w[i];
  }
}

/**
 * Editor for the Matrix field of the config mode
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_mat() {
  if (key == CODED) {
    if (keyCode == UP) {
      if (config_mat_cind[1] == 0) {
        config_editor_field = config_field.DIM;
        config_editor_entry = config_dim_str.length();
        update_mat();
      } else {
        config_mat_cind[1]--;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = min(config_editor_entry, s.length());
      }
    } else if (keyCode == DOWN) {
      if (config_mat_cind[1] == state_dim - 1) {
        config_editor_field = config_field.VEC;
        config_vec_cind     = 0;
        config_editor_entry = config_vec_str[config_vec_cind].length();
        update_mat();
      } else {
        config_mat_cind[1]++;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = min(config_editor_entry, s.length());
      }
    } else if (keyCode == LEFT) {
      if (config_editor_entry == 0
      &&  config_mat_cind[0]  == 0
      &&  config_mat_cind[1]  == 0) {
        config_editor_field = config_field.DIM;
        config_editor_entry = config_dim_str.length();
        update_mat();
      } else if (config_editor_entry == 0
             &&  config_mat_cind[0]  == 0) {
        config_mat_cind[0] = state_dim - 1;
        config_mat_cind[1]--;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
      } else if (config_editor_entry == 0) {
        config_mat_cind[0]--;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
      } else {
        config_editor_entry--;
      }
    } else if (keyCode == RIGHT) {
      int end = config_mat_str[config_mat_cind[0]][config_mat_cind[1]].length();
      if (config_editor_entry == end
      &&  config_mat_cind[0]  == state_dim - 1
      &&  config_mat_cind[1]  == state_dim - 1) {
        config_editor_field = config_field.VEC;
        config_editor_entry = 0;
        config_vec_cind     = 0;
        update_mat();
      } else if (config_editor_entry == end
             &&  config_mat_cind[0]  == state_dim - 1) {
        config_mat_cind[0]  = 0;
        config_mat_cind[1]++;
        config_editor_entry = 0;
      } else if (config_editor_entry == end) {
        config_mat_cind[0]++;
        config_editor_entry = 0;
      } else {
        config_editor_entry++;
      }
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      if (config_mat_cind[0] == state_dim - 1
      &&  config_mat_cind[1] == state_dim - 1) {
        config_editor_field = config_field.VEC;
        config_vec_cind     = 0;
        config_editor_entry = config_vec_str[config_vec_cind].length();
        update_mat();
      } else if (config_mat_cind[0] == state_dim - 1) {
        config_mat_cind[0] = 0;
        config_mat_cind[1]++;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
      } else {
        config_mat_cind[0]++;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
      }
    } else if (key == 8) {
      // Backspace
      String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
      if (config_editor_entry > 0) {
        s = string_remove(s, config_editor_entry - 1);
        config_mat_str[config_mat_cind[0]][config_mat_cind[1]] = s;
        config_editor_entry--;
        update_strmat();
      }
    } else if (key == 127) {
      // Delete
      String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
      if (config_editor_entry < s.length()) {
        s = string_remove(s, config_editor_entry);
        config_mat_str[config_mat_cind[0]][config_mat_cind[1]] = s;
        update_strmat();
      }
    } else if ('0' <= key && key <= '9' || key == '.') {
      String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
      if (key == '.' && s.indexOf('.') > -1)
        return;
      s = string_append(s, config_editor_entry, key);
      config_mat_str[config_mat_cind[0]][config_mat_cind[1]] = s;
      config_editor_entry++;
      update_strmat();
    }
  }
}

/**
 * Updates state_vec based on user input
 *
 * Called after the user leaves the Vector field.
 */
void update_vec() {
  state_vec = new float[1][state_dim];
  for (int i = 0; i < state_dim; i++) {
    if (config_vec_str[i].equals("") || config_vec_str[i].equals("."))
      config_vec_str[i] = "0";
    state_vec[0][i] = float(config_vec_str[i]);
    if (state_vec[0][i] == int(state_vec[0][i]))
      config_vec_str[i] = str(int(state_vec[0][i]));
    else
      config_vec_str[i] = str(state_vec[0][i]);
  }
} 

/**
 * Editor for the Vector field of the config mode
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_vec() {
  if (key == CODED) {
    if (keyCode == UP) {
      if (config_vec_cind == 0) {
        config_editor_field = config_field.MAT;
        config_mat_cind[0]  = 0;
        config_mat_cind[1]  = state_dim - 1;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
        update_vec();
      } else {
        config_vec_cind--;
        String s = config_vec_str[config_vec_cind];
        config_editor_entry = min(config_editor_entry, s.length());
      }
    } else if (keyCode == DOWN) {
      if (config_vec_cind == state_dim - 1) {
        config_editor_field = config_field.DT;
        config_editor_entry = config_dt_str.length();
        update_vec();
      } else {
        config_vec_cind++;
        String s = config_vec_str[config_vec_cind];
        config_editor_entry = min(config_editor_entry, s.length());
      }
    } else if (keyCode == LEFT) {
      if (config_editor_entry == 0 && config_vec_cind == 0) {
        config_editor_field = config_field.MAT;
        config_mat_cind[0]  = state_dim - 1;
        config_mat_cind[1]  = state_dim - 1;
        String s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
        config_editor_entry = s.length();
        update_vec();
      } else if (config_editor_entry == 0) {
        config_vec_cind--;
        config_editor_entry = config_vec_str[config_vec_cind].length();
      } else {
        config_editor_entry--;
      }
    } else if (keyCode == RIGHT) {
      int end = config_vec_str[config_vec_cind].length();
      if (config_editor_entry == end && config_vec_cind == state_dim - 1) {
        config_editor_field = config_field.DT;
        config_editor_entry = 0;
        update_vec();
      } else if (config_editor_entry == end) {
        config_vec_cind++;
        config_editor_entry = 0;
      } else {
        config_editor_entry++;
      }
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      if (config_vec_cind == state_dim - 1) {
        config_editor_field = config_field.DT;
        config_editor_entry = config_dt_str.length();
        update_vec();
      } else {
        config_vec_cind++;
        config_editor_entry = config_vec_str[config_vec_cind].length();
      }
    } else if (key == 8) {
      // Backspace
      String s = config_vec_str[config_vec_cind];
      if (config_editor_entry > 0) {
        s = string_remove(s, config_editor_entry - 1);
        config_vec_str[config_vec_cind] = s;
        config_editor_entry--;
      }
    } else if (key == 127) {
      // Delete
      String s = config_vec_str[config_vec_cind];
      if (config_editor_entry < s.length()) {
        s = string_remove(s, config_editor_entry);
        config_vec_str[config_vec_cind] = s;
      }
    } else if ('0' <= key && key <= '9' || key == '.') {
      String s = config_vec_str[config_vec_cind];
      if (key == '.' && s.indexOf('.') > -1)
        return;
      s = string_append(s, config_editor_entry, key);
      config_vec_str[config_vec_cind] = s;
      config_editor_entry++;
    }
  }
}

/**
 * Updates state_dt based on user input
 *
 * Called after the user leaves the Time Step Size field. There is a default
 * value of 0.025 and a max value of 1.
 */
void update_dt() {
  if (config_dt_str.equals("") || config_dt_str.equals("."))
    config_dt_str = "0.025";
  state_dt = float(config_dt_str);
  if (state_dt >= 1) {
    config_dt_str = "1";
    state_dt = 1;
  } else {
    config_dt_str = str(state_dt);
  }
}

/**
 * Editor for the Time Step Size field
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_dt() {
  if (key == CODED) {
    if (keyCode == UP) {
      config_editor_field = config_field.VEC;
      config_vec_cind     = state_dim - 1;
      config_editor_entry = config_vec_str[config_vec_cind].length();
      update_dt();
    } else if (keyCode == DOWN) {
      config_editor_field = config_field.NTOTAL;
      config_editor_entry = config_ntotal_str.length();
      update_dt();
    } else if (keyCode == LEFT) {
        if (config_editor_entry == 0) {
          config_editor_field = config_field.VEC;
          config_vec_cind     = state_dim - 1;
          config_editor_entry = config_vec_str[config_vec_cind].length();
          update_dt();
        } else
          config_editor_entry--;
    } else if (keyCode == RIGHT) {
        if (config_editor_entry == config_dt_str.length()) {
          config_editor_field = config_field.NTOTAL;
          config_editor_entry = 0;
          update_dt();
        } else
          config_editor_entry++;
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      config_editor_field = config_field.NTOTAL;
      config_editor_entry = config_ntotal_str.length();
      update_dt();
    } else if (key == 8 && config_editor_entry > 0) {
      // Backspace
      config_dt_str = string_remove(config_dt_str, config_editor_entry - 1);
      config_editor_entry--;
    } else if (key == 127 && config_editor_entry < config_dt_str.length()) {
      // Delete
      config_dt_str = string_remove(config_dt_str, config_editor_entry);
    } else if ('0' <= key && key <= '9' || key == '.') {
      if (key == '.' && config_dt_str.indexOf('.') > -1)
        return;
      config_dt_str = string_append(config_dt_str, config_editor_entry, key);
      config_editor_entry++;
    }
  }
}

/**
 * Updates state_ntotal based on user input
 *
 * Called after the user leaves the Total Nodes field. There is a default
 * value of 50.
 */
void update_ntotal() {
  if (config_ntotal_str.equals(""))
    config_ntotal_str = "50";
  state_ntotal = int(config_ntotal_str);
  config_ntotal_str = str(state_ntotal);
}

/**
 * Editor for the Total Nodes field
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_ntotal() {
  if (key == CODED) {
    if (keyCode == UP) {
      config_editor_field = config_field.DT;
      config_editor_entry = config_dt_str.length();
      update_ntotal();
    } else if (keyCode == DOWN) {
      config_editor_field = config_field.NSPREAD;
      config_editor_entry = config_nspread_str.length();
      update_ntotal();
    } else if (keyCode == LEFT) {
        if (config_editor_entry == 0) {
          config_editor_field = config_field.DT;
          config_editor_entry = config_dt_str.length();
          update_ntotal();
        } else
          config_editor_entry--;
    } else if (keyCode == RIGHT) {
        if (config_editor_entry == config_ntotal_str.length()) {
          config_editor_field = config_field.NSPREAD;
          config_editor_entry = 0;
          update_ntotal();
        } else
          config_editor_entry++;
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      config_editor_field = config_field.NSPREAD;
      config_editor_entry = config_nspread_str.length();
      update_ntotal();
    } else if (key == 8) {
      // Backspace
      String s = config_ntotal_str;
      if (config_editor_entry > 0) {
        config_ntotal_str = string_remove(s, config_editor_entry - 1);
        config_editor_entry--;
      }
    } else if (key == 127) {
      // Delete
      String s = config_ntotal_str;
      if (config_editor_entry < s.length())
        config_ntotal_str = string_remove(s, config_editor_entry);
    } else if ('0' <= key && key <= '9') {
      String s = config_ntotal_str;
      config_ntotal_str = string_append(s, config_editor_entry, key);
      config_editor_entry++;
    }
  }
}

/**
 * Updates state_nspread based on user input
 *
 * Called after the user leaves the Node Spread Ratio field. There is a default
 * value of 0.5.
 */
void update_nspread() {
  if (config_nspread_str.equals("") || config_nspread_str.equals("."))
    config_nspread_str = "0.5";
  state_nspread = float(config_nspread_str);
  if (state_nspread == int(state_nspread))
    config_nspread_str = str(int(state_nspread));
  else
    config_nspread_str = str(state_nspread);
}

/**
 * Editor for the Node Spread Ratio field
 *
 * Changes the text displayed to the user but does not change any underlying
 * variables until exitting.
 */
void edit_nspread() {
  if (key == CODED) {
    if (keyCode == UP) {
      config_editor_field = config_field.NTOTAL;
      config_editor_entry = config_ntotal_str.length();
      update_nspread();
    } else if (keyCode == DOWN) {
      config_editor_field = config_field.CONT;
      config_editor_entry = 0;
      update_nspread();
    } else if (keyCode == LEFT) {
        if (config_editor_entry == 0) {
          config_editor_field = config_field.NTOTAL;
          config_editor_entry = config_ntotal_str.length();
          update_nspread();
        } else
          config_editor_entry--;
    } else if (keyCode == RIGHT) {
        if (config_editor_entry == config_nspread_str.length()) {
          config_editor_field = config_field.CONT;
          config_editor_entry = 0;
          update_nspread();
        } else
          config_editor_entry++;
    }
  } else {
    if (key == 9 || key == 10) {
      // Tab or Enter
      config_editor_field = config_field.CONT;
      config_editor_entry = 0;
      update_nspread();
    } else if (key == 8) {
      // Backspace
      String s = config_nspread_str;
      if (config_editor_entry > 0) {
        config_nspread_str = string_remove(s, config_editor_entry - 1);
        config_editor_entry--;
      }
    } else if (key == 127) {
      // Delete
      String s = config_nspread_str;
      if (config_editor_entry < s.length())
        config_nspread_str = string_remove(s, config_editor_entry);
    } else if ('0' <= key && key <= '9' || key == '.') {
      String s = config_nspread_str;
      if (key == '.' && s.indexOf('.') > -1)
        return;
      config_nspread_str = string_append(s, config_editor_entry, key);
      config_editor_entry++;
    }
  }
}

/**
 * Editor for the Continue field
 */
void edit_cont() {
  if (key == CODED) {
    if (keyCode == UP || keyCode == LEFT) {
      config_editor_field = config_field.NSPREAD;
      config_editor_entry = config_nspread_str.length();
    } else if (keyCode == DOWN || keyCode == RIGHT) {
      config_editor_field = config_field.DIM;
      config_editor_entry = 0;
    }
  } else {
    if (key == 9) {
      // Tab
      config_editor_field = config_field.DIM;
      config_editor_entry = config_dim_str.length();
    } else if (key == 10) {
      // Enter
      mode = app_mode.MAIN;
      config_editor_field = config_field.DIM;
      config_editor_entry = 0;
      config_scroll_x     = 0;
      config_scroll_y     = 0;
      setup_main();
    }
  }
}

/**
 * Updates the scroll values by delegating based on the current mode
 */
void update_scroll() {
  float xi, xf, y;
  String s;
  float w = width  - text_margin - 3 * text_size;
  float h = height - text_margin - 3 * text_size;
  switch (config_editor_field) {
  case DIM:
    xi = config_dim_x;
    xf = config_dim_x;
    y  = config_dim_y;
    s  = config_dim_str;
    break;
  case MAT:
    xi = config_mat_x;
    xf = config_mat_coords[config_mat_cind[0]][config_mat_cind[1]][0];
    y  = config_mat_coords[config_mat_cind[0]][config_mat_cind[1]][1];
    s  = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
    break;
  case VEC:
    xi = config_vec_x;
    xf = config_vec_x;
    y  = config_vec_coords[config_vec_cind];
    s  = config_vec_str[config_vec_cind];
    break;
  case DT:
    xi = config_dt_x;
    xf = config_dt_x;
    y  = config_dt_y;
    s  = config_dt_str;
    break;
  case NTOTAL:
    xi = config_ntotal_x;
    xf = config_ntotal_x;
    y  = config_ntotal_y;
    s  = config_ntotal_str;
    break;
  case NSPREAD:
    xi = config_nspread_x;
    xf = config_nspread_x;
    y  = config_nspread_y;
    s  = config_nspread_str;
    break;
  case CONT:
    config_scroll_x = 0;
    config_scroll_y = max(config_cont_y - (height - text_margin - text_size), 0);
    return;
  default:
    throw new RuntimeException("Invalid field");
  }
  s  = s.substring(0, config_editor_entry);
  xf += textWidth(s);
  if (xf - config_scroll_x < xi)
    config_scroll_x = xf - xi;
  else if (xf - config_scroll_x > w)
    config_scroll_x = xf - w;
  if (config_editor_field == config_field.DIM)
    config_scroll_y = 0;
  else if (y - config_scroll_y < config_mat_y)
    config_scroll_y = y - config_mat_y;
  else if (y - config_scroll_y > h)
    config_scroll_y = y - h;  
}

/**
 * Handle key presses by delegating based on the current field being editted
 */
void keyPressed_config() {
  switch (config_editor_field) {
  case DIM:
    edit_dim();
    break;
  case MAT:
    edit_mat();
    break;
  case VEC:
    edit_vec();
    break;
  case DT:
    edit_dt();
    break;
  case NTOTAL:
    edit_ntotal();
    break;
  case NSPREAD:
    edit_nspread();
    break;
  case CONT:
    edit_cont();
    break;
  default:
    throw new RuntimeException("Invalid field");
  }
  update_scroll();
}

/**
 * Helper function for handling jumping during mouse presses
 *
 * Returns false if no jump occured.
 * Returns true  if a  jump occured.
 *
 */
boolean mousePressed_config_helper() {
  mouseX += config_scroll_x;
  mouseY += config_scroll_y;
  if (mouseY < config_dim_y || mouseY > config_cont_y + text_size)
    return false;
  else if (mouseY < config_mat_y) {
    if (inbtw(config_dim_x, config_dim_x + textWidth(config_dim_str), mouseX)
    &&  mouseY <= config_dim_y + text_size) {
      config_editor_field = config_field.DIM;
      config_editor_entry = config_dim_str.length();
    } else return false;
  } else if (mouseY < config_vec_y) {
    int i = -1;
    for (int k = 0; k < state_dim; k++) {
      float xstart = config_mat_coords[k][0][0];
      if (inbtw(xstart, xstart + config_mat_w[k], mouseX)) {
        i = k;
        break;
      }
    }
    if (i == -1) return false;
    int j = int((mouseY - config_mat_y) / (2 * text_size));
    if (0 <= j && j < state_dim
    &&  mouseY <= config_mat_coords[i][j][1] + text_size) {
      config_editor_field = config_field.MAT;
      config_mat_cind[0]  = i;
      config_mat_cind[1]  = j;
      config_editor_entry = config_mat_str[i][j].length();
    } else return false;
  } else if (mouseY < config_dt_y) {
    int i = int((mouseY - config_vec_y) / (2 * text_size));
    float w = max(textWidth("_"), textWidth(config_vec_str[i]));
    if (inbtw(config_vec_x, config_vec_x + w, mouseX)
    &&  mouseY <= config_vec_coords[i] + text_size) {
      config_editor_field = config_field.VEC;
      config_vec_cind     = i;
      config_editor_entry = config_vec_str[i].length();
    } else return false;
  } else if (mouseY < config_ntotal_y) {
    float w = max(textWidth("_"), textWidth(config_dt_str));
    if (inbtw(config_dt_x, config_dt_x + w, mouseX)
    &&  mouseY < config_dt_y + text_size) {
      config_editor_field = config_field.DT;
      config_editor_entry = config_dt_str.length();
    } else return false;
  } else if (mouseY < config_nspread_y) {
    float w = max(textWidth("_"), textWidth(config_ntotal_str));
    if (inbtw(config_ntotal_x, config_ntotal_x + w, mouseX)
    &&  mouseY < config_ntotal_y + text_size) {
      config_editor_field = config_field.NTOTAL;
      config_editor_entry = config_ntotal_str.length();
    } else return false;
  } else if (mouseY < config_cont_y) {
    float w = max(textWidth("_"), textWidth(config_nspread_str));
    if (inbtw(config_nspread_x, config_nspread_x + w, mouseX)
    &&  mouseY < config_nspread_y + text_size) {
      config_editor_field = config_field.NSPREAD;
      config_editor_entry = config_nspread_str.length();
    } else return false;
  } else {
    config_editor_field = config_field.CONT;
    config_editor_entry = 0;
  }
  return true;
}

/**
 * If the mouse is pressed at an entry field then jump to that field and update
 * the old field
 *
 * If the user has changed the state_dim field then just update the dimension
 * and all dependent variables.
 */
void mousePressed_config() {
  if (config_editor_field == config_field.DIM
  &&  state_dim != int(config_dim_str)) {
    update_dim();
    return;
  }
  config_field prev = config_editor_field;
  if (!mousePressed_config_helper())
    return;

  if (prev != config_editor_field) {
    switch (prev) {
    case DIM:
      update_dim();
      break;
    case MAT:
      update_mat();
      break;
    case VEC:
      update_vec();
      break;
    case DT:
      update_dt();
      break;
    case NTOTAL:
      update_ntotal();
      break;
    case NSPREAD:
      update_nspread();
      break;
    default:
      throw new RuntimeException("Invalid field");
    }
  }

  if (config_editor_field != config_field.CONT)
    update_scroll();
  else {
    mode = app_mode.MAIN;
    config_editor_field = config_field.DIM;
    config_editor_entry = 0;
    config_scroll_x = 0;
    config_scroll_y = 0;
    update_vec();
    setup_main();
  }
}

/**
 * If the mouse wheel is moved then vertically scroll accordingly
 *
 * Note: The horizontal scroll is unchangable as Processing cannot detect
 * horizontal wheel movement.
 */
void mouseWheel_config(MouseEvent e) {
  config_scroll_y += height / 50 * e.getCount();
  if (text_margin - config_scroll_y > text_margin)
    config_scroll_y = 0;
  else if (config_cont_y + text_size - config_scroll_y < height - text_margin)
    config_scroll_y = config_cont_y + text_size - (height - text_margin);
}


/********** Event Functions for Main **********/

/**
 * If rewinding an iteration then pull from history
 */
void anim_back() {
  if (iters == 0)
    return;
  iters--;
  state_vec = hist_states.get(iters);
}

/**
 * If forwarding an iteration then pull from history if possible; otherwise
 * calculate the number of nodes traveling along each edge
 */
void anim_forward() {
  if (iters < hist_ninit.size()) {
    node_init   = hist_ninit.get(iters);
    node_finl   = hist_nfinl.get(iters);
    node_counts = hist_ncnts.get(iters);
    node_len    = hist_nlen.get(iters);
    state_vec   = hist_states.get(iters + 1);
  } else {
    float[][] v = matrixMultiply(state_mat, state_vec);
    int n       = state_dim * state_dim - state_dim;
    node_init   = new float[n][2];
    node_finl   = new float[n][2];
    node_counts = new float[n];
    node_len    = 0;
    float ratio = weightRatio(state_vec, state_ntotal);
    for (int i = 0; i < state_dim; i++) {
      for (int j = 0; j < state_dim; j++) {
        if (i != j && state_mat[i][j] != 0 && state_vec[0][i] != 0) {
          node_init[node_len]   = vert_pos[i];
          node_finl[node_len]   = vert_pos[j];
          node_counts[node_len] = state_mat[i][j] * state_vec[0][i] * ratio;
          node_len++;
        }
      }
    }
    hist_ninit.add(node_init);
    hist_nfinl.add(node_finl);
    hist_ncnts.add(node_counts);
    hist_nlen.append(node_len);
    hist_states.add(v);
    state_vec = v;
  }
  iters++;
  isAnimating = true;
}

/**
 * Key press handler
 *
 * If the RIGHT key is pressed then increment the iteration count.
 * If the LEFT  key is pressed then decrement the iteration count.
 */
void keyPressed_main() {
  if (isAnimating) {
    if (key == CODED && keyCode == RIGHT) {
      t = 0;
      isAnimating = false;
    }
  } else if (key == CODED) {
    if (keyCode == RIGHT)
      anim_forward();
    else if (keyCode == LEFT)
      anim_back();
  }
}

/**
 * Mouse press handler
 *
 * If MB1/LEFT  is pressed then increment the iteration count.
 * If MB2/RIGHT is pressed then decrement the iteration count.
 */
void mousePressed_main() {
  if (isAnimating) {
    if (mouseButton == LEFT) {
      t = 0;
      isAnimating = false;
    }
  } else {
    if (mouseButton == LEFT)
      anim_forward();
    else if (mouseButton == RIGHT)
      anim_back();
  }
}

/**
 * If the mouse wheel is moved then change iteration accordingly
 *
 * Upward   movement is reflected by incrementing the iteration count.
 * Downward movement is reflected by decrementing the iteration count.
 */
void mouseWheel_main(MouseEvent e) {
  if (isAnimating) {
    if (e.getCount() > 0) {
      t = 0;
      isAnimating = false;
    }
  } else {
    if (e.getCount() > 0)
      anim_forward();
    else if (e.getCount() < 0)
      anim_back();
  }
}


/********** Event Functions **********/

/** Key Press Handler */
void keyPressed() {
  if (key == 'm' || key == 'M')
    mode = app_mode.MENU;
  else if (key == 'h' || key == 'H')
    mode = app_mode.HELP;
  else if (key == 'c' || key == 'C')
    mode = app_mode.CONFIG;
  else if (mode == app_mode.CONFIG)
    keyPressed_config();
  else if (mode == app_mode.MAIN)
    keyPressed_main();
}

/** Mouse Press Handler */
void mousePressed() {
  if (mode == app_mode.CONFIG)
    mousePressed_config();
  else if (mode == app_mode.MAIN)
    mousePressed_main();
}

/** Mouse Wheel Handler */
void mouseWheel(MouseEvent e) {
  if (mode == app_mode.CONFIG)
    mouseWheel_config(e);
  else if (mode == app_mode.MAIN)
    mouseWheel_main(e);
}


/********** Draw Functions **********/

/**
 * Draws the menu mode
 */
void draw_menu() {
  background(0);
  textAlign(CENTER, BOTTOM);
  fill(255);
  textSize(text_size * 3);
  text("Welcome!", width / 2, height / 2);
  textSize(text_size);
  text("Press C for Configurations", width / 2, height / 2 + 5 * text_size);
  text("Press H for Help", width / 2, height / 2 + 6 * text_size);
}

/**
 * Draws the help mode
 */
void draw_help() {
  background(0);
  textAlign(LEFT, TOP);
  fill(0, 255, 0);
  text(help_hdr, text_margin, text_margin);
  float x = text_indent;
  float y = text_margin + 2 * text_size;
  float w = (width  - text_margin) - x;
  float h = (height - text_margin) - y;
  text(help_msg, x, y, w, h);
}

/**
 * Draws the current field/entry indicator in config mode
 */
void draw_config_curr() {
  stroke(255, 0, 0);
  fill(255, 0, 0);
  float x, y;
  String s;
  if (config_editor_field != config_field.CONT) {
    switch (config_editor_field) {
    case DIM:
      x = config_dim_x;
      y = config_dim_y;
      s = config_dim_str;
      break;
    case MAT:
      x = config_mat_coords[config_mat_cind[0]][config_mat_cind[1]][0];
      y = config_mat_coords[config_mat_cind[0]][config_mat_cind[1]][1];
      s = config_mat_str[config_mat_cind[0]][config_mat_cind[1]];
      break;
    case VEC:
      x = config_vec_x;
      y = config_vec_coords[config_vec_cind];
      s = config_vec_str[config_vec_cind];
      break;
    case DT:
      x = config_dt_x;
      y = config_dt_y;
      s = config_dt_str;
      break;
    case NTOTAL:
      x = config_ntotal_x;
      y = config_ntotal_y;
      s = config_ntotal_str;
      break;
    case NSPREAD:
      x = config_nspread_x;
      y = config_nspread_y;
      s = config_nspread_str;
      break;
    default:
      throw new RuntimeException("Invalid field");
    }
    s = s.substring(0, config_editor_entry);
    x += textWidth(s);
    x -= config_scroll_x;
    y -= config_scroll_y;
    line(x, y, x, y + text_size);
  } else {
    x = text_margin   - config_scroll_x;
    y = config_cont_y - config_scroll_y;
    s = "Press here to continue.";
    rect(x, y, textWidth(s), text_size);
  }
}

/**
 * Helper function to draw one-line fields in the config mode
 */
void draw_config_helper(String title, String s, float x, float y) {
  text(title, text_indent - config_scroll_x, y - config_scroll_y);
  if (s.equals(""))
    text("_", x - config_scroll_x, y - config_scroll_y);
  else
   text(s, x - config_scroll_x, y - config_scroll_y);
}

/**
 * Draws the config mode
 */
void draw_config() {
  float sx = config_scroll_x;
  float sy = config_scroll_y;
  background(0);
  textAlign(LEFT, TOP);
  fill(0, 0, 255);
  text("Please enter the parameters.", text_margin - sx, text_margin - sy);
  draw_config_helper("Dimension:", config_dim_str, config_dim_x, config_dim_y);
  text("Matrix:", text_indent - sx, config_mat_y - sy);
  for (int i = 0; i < state_dim; i++) {
    for (int j = 0; j < state_dim; j++) {
      String s;
      if (config_mat_str[i][j].equals(""))
        s = "_";
      else
        s = config_mat_str[i][j];
      text(s, config_mat_coords[i][j][0] - sx, config_mat_coords[i][j][1] - sy);
    }
  }
  text("Initial State:", text_indent - sx, config_vec_y - sy);
  for (int i = 0; i < state_dim; i++) {
    if (config_vec_str[i].equals(""))
      text("_", config_vec_x - sx, config_vec_coords[i] - sy);
    else
      text(config_vec_str[i], config_vec_x - sx, config_vec_coords[i] - sy);
  }
  draw_config_helper("Time Step Size:",
                     config_dt_str,
                     config_dt_x,
                     config_dt_y);
  draw_config_helper("Total Nodes:",
                     config_ntotal_str,
                     config_ntotal_x,
                     config_ntotal_y);
  draw_config_helper("Node Spread Ratio:",
                     config_nspread_str,
                     config_nspread_x,
                     config_nspread_y);
  draw_config_curr();
  fill(0, 0, 255);
  text("Press here to continue.", text_margin - sx, config_cont_y - sy);
}

/**
 * Draws the first set of state information during main mode
 *
 * Draws the edges "underneath" everything else.
 */
void draw_info_under() {
  background(0);
  stroke(255);
  for (int i = 0; i < state_dim; i++) {
    for (int j = 0; j < state_dim; j++) {
      if (state_mat[i][j] != 0)
        line(vert_pos[i][0], vert_pos[i][1], vert_pos[j][0], vert_pos[j][1]);
    }
  }
}

/**
 * Draws the current frame of the animation during main mode
 */
void draw_anim() {
  textAlign(RIGHT, BOTTOM);
  fill(0, 255, 0);
  textSize(text_size);
  text("Loading...", width, height);
  stroke(255);
  if (node_len == 0) {
    isAnimating = false;
    return;
  }
  for (int i = 0; i < node_len; i++) {
    float dx = node_finl[i][0] - node_init[i][0];
    float dy = node_finl[i][1] - node_init[i][1];
    float cx = node_init[i][0] + t * dx;
    float cy = node_init[i][1] + t * dy;
    fill(vert_colors[i][0], vert_colors[i][1], vert_colors[i][2]);
    for (int j = 0; j < node_counts[i]; j++) {
      float n = state_nspread * j / (node_counts[i] - 1);
      float x = cx - n * dx;
      float y = cy - n * dy;
      if (inbtw(node_init[i][0], node_finl[i][0], x)
      &&  inbtw(node_init[i][1], node_finl[i][1], y))
        ellipse(x, y, node_d, node_d);
    }
  }
  t += state_dt;
  if (t > 1 + state_nspread) {
    t = 0;
    isAnimating = false;
  }
}

/**
 * Draws the last set of state information (everything else) during main mode
 *
 * Draws the rest of the information "above" the edges and animation by being
 * called last.
 */
void draw_info_above() {
  noStroke();
  textAlign(CENTER, CENTER);
  textSize(vert_textsize);
  for (int i = 0; i < state_dim; i++) {
    fill(255);
    ellipse(vert_pos[i][0], vert_pos[i][1], vert_d, vert_d);
    fill(0, 0, 255);
    text(str(i + 1), vert_pos[i][0], vert_pos[i][1]);
  }
  fill(255, 0, 0);
  for (int i = 0; i < state_dim; i++) {
    String s = str(state_vec[0][i]);
    text(s, vert_tpos[i][0], vert_tpos[i][1]);
  }
  textAlign(LEFT, BOTTOM);
  fill(0, 255, 0);
  textSize(text_size);
  text("Iterations: " + str(iters), 0, height);
}

/**
 * Draws the main mode
 */
void draw_main() {
  draw_info_under();
  if (isAnimating)
    draw_anim();
  draw_info_above();
}

/**
 * Draw by delegating according to mode
 */
void draw() {
  switch (mode) {
  case MENU:
    draw_menu();
    break;
  case HELP:
    draw_help();
    break;
  case CONFIG:
    draw_config();
    break;
  case MAIN:
    draw_main();
    break;
  default:
    throw new RuntimeException("Invalid mode");
  }
}