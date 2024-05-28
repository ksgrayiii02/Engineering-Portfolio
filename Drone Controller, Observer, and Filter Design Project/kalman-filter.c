// -------------------- KF Filter --------------------
void predict()
{
  float temp1[12][12] = {}; // temporary matrix
  float temp2[12][12] = {}; // temporary matrix

  // x = Fx + Gu
  matMultiply(12,12,12,1,F,x,temp1);
  matMultiply(12,4,4,1,G,u,temp2);
  matAdd(12,1,temp1,temp2,x);
  
  // P = FPF.T + Q
  matMultiply(12,12,12,12,F,P,temp1);
  matMultiply(12,12,12,12,temp1,FT,temp2);
  matAdd(12,12,temp2,Q,P);      

}

void update()
{
  float K[12][3] = {};
  float KT[3][12] = {};

  float temp1[12][3] = {}; // temporary matrix
  float temp2[3][12] = {}; // temporary matrix
  float temp3[3][3] = {};  // temporary matrix
  float temp4[12][12] = {};
  float temp5[12][12] = {};
  float temp6[12][12] = {};

  float R[3][3] = {
    {10.969344f,0.0f,0.0f},
    {0.0f,27.888961f,0.0f},
    {0.0f,0.0f,9e-06f},
  };

  float I[12][12] = {
    {1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f},
    {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f}
  };

  float dh_dx[3][12] = {
    {0.0f,0.0f,-4.09255567950588f*w_y/o_z - (-4.09255567950588f*o_z*w_y + 4.09255567950588f*v_x)/powf(o_z,2),0.0f,0.0f,0.0f,4.09255567950588f/o_z,0.0f,0.0f,0.0f,-4.09255567950588f,0.0f},
    {0.0f,0.0f,4.09255567950588f*w_x/o_z - (4.09255567950588f*o_z*w_x + 4.09255567950588f*v_y)/powf(o_z,2),0.0f,0.0f,0.0f,0.0f,4.09255567950588f/o_z,0.0f,4.09255567950588f,0.0f,0.0f},
    {0.0f,0.0f,1.0f/(cosf(phi)*cosf(theta)),0.0f,o_z*sinf(theta)/(cosf(phi)*powf(cosf(theta),2)),o_z*sinf(phi)/(powf(cosf(phi),2)*cosf(theta)),0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
  };

  float dh_dx_T[12][3] = {};
  transpose(3,12,dh_dx,dh_dx_T);

  // K = PH.T*(HP*H.T+R).inv
  matMultiply(12,12,12,3,P,dh_dx_T,temp1);

  matMultiply(3,12,12,12,dh_dx,P,temp2);
  matMultiply(3,12,12,3,temp2,dh_dx_T,temp3);
  matAdd(3,3,temp3,R,temp3);

  inverse3x3(temp3);

  matMultiply(12,3,3,3,temp1,temp3,K);

  // x_new = x_old + K
  matAdd(12,3,x,K,x);

  // P_new = (I-KH)P*(I-KH).T+KRK.T
  matMultiply(12,3,3,12,K,dh_dx,temp4);
  matSubtract(12,12,I,temp4,temp4); // temp4 = (I-KH)
  
  transpose(12,12,temp4,temp5); // temp5 = temp4.T = (I-KH).T

  matMultiply(12,12,12,12,temp4,P,temp6);  // temp6 = (I-KH)P
  matMultiply(12,12,12,12,temp6,temp5,temp4); // temp4 = temp6 * temp5 = (I-KH)P * (I-KH).T

  matMultiply(12,3,3,3,K,R,temp1);
  transpose(12,3,K,KT);
  matMultiply(12,3,3,12,temp1,KT,temp5); // temp5 = KRK.T

  matAdd(12,12,temp4,temp5,P);

}
// -------------------- KF Filter --------------------

void controllerAE483(control_t *control,
                     setpoint_t *setpoint,
                     const sensorData_t *sensors,
                     const state_t *state,
                     const uint32_t tick)
{
  if (RATE_DO_EXECUTE(ATTITUDE_RATE, tick)) {
    // Everything in here runs at 500 Hz

    // Desired position
    o_x_des = setpoint->position.x;
    o_y_des = setpoint->position.y;
    o_z_des = setpoint->position.z;

    // Measurements
    w_x = radians(sensors->gyro.x);
    w_y = radians(sensors->gyro.y);
    w_z = radians(sensors->gyro.z);
    a_z = g * sensors->acc.z;
    n_x = flow_dpixelx;
    n_y = flow_dpixely;
    r = tof_distance;

    if (reset_observer) {
      o_x = 0.0f;
      o_y = 0.0f;
      o_z = 0.0f;
      psi = 0.0f;
      theta = 0.0f;
      phi = 0.0f;
      v_x = 0.0f;
      v_y = 0.0f;
      v_z = 0.0f;
      reset_observer = false;
    }

    // State estimates
    if (use_observer) {
    
      // Compute each element of:
      // 
      //   C x + D u - y
      // 
      // FIXME: your code goes here
      n_x_err = k_flow * ((v_x / o_z_eq) - w_y) - n_x;  // <-- FIXME
      n_y_err = k_flow * ((v_y / o_z_eq) + w_x) - n_y;  // <-- FIXME
      r_err = o_z - r;    // <-- FIXME
      
      // Update estimates
      // FIXME: your code goes here
      o_x += dt * (v_x);
      o_y += dt * (v_y); 
      o_z += dt * (v_z - 19.922768f*r_err);
      psi += dt * (w_z);
      theta += dt * (w_y - 0.004591f*n_x_err);
      phi += dt * (w_x - -0.031130f*n_y_err);
      v_x += dt * (g*theta - 0.108703f*n_x_err); 
      v_y += dt * (-g*phi - 0.273845f*n_y_err); 
      v_z += dt * (a_z-g - 93.333333f*r_err); 

      x[0][0] = o_x; // UPDATE FOR KF
      x[1][0] = o_y;
      x[2][0] = o_z;
      x[3][0] = psi;
      x[4][0] = theta;
      x[5][0] = phi;
      x[6][0] = v_x;
      x[7][0] = v_y;
      x[8][0] = v_z;
      
    } else {
      o_x = state->position.x;
      o_y = state->position.y;
      o_z = state->position.z;
      psi = radians(state->attitude.yaw);
      theta = - radians(state->attitude.pitch);
      phi = radians(state->attitude.roll);
      v_x = state->velocity.x*cosf(psi)*cosf(theta) + state->velocity.y*sinf(psi)*cosf(theta) - state->velocity.z*sinf(theta);
      v_y = state->velocity.x*(sinf(phi)*sinf(theta)*cosf(psi) - sinf(psi)*cosf(phi)) + state->velocity.y*(sinf(phi)*sinf(psi)*sinf(theta) + cosf(phi)*cosf(psi)) + state->velocity.z*sinf(phi)*cosf(theta);
      v_z = state->velocity.x*(sinf(phi)*sinf(psi) + sinf(theta)*cosf(phi)*cosf(psi)) + state->velocity.y*(-sinf(phi)*cosf(psi) + sinf(psi)*sinf(theta)*cosf(phi)) + state->velocity.z*cosf(phi)*cosf(theta);

      x[0][0] = o_x; // UPDATE FOR KF
      x[1][0] = o_y;
      x[2][0] = o_z;
      x[3][0] = psi;
      x[4][0] = theta;
      x[5][0] = phi;
      x[6][0] = v_x;
      x[7][0] = v_y;
      x[8][0] = v_z;
    }

    if (setpoint->mode.z == modeDisable) {
      // If there is no desired position, then all
      // motor power commands should be zero

      u[0][0] = 0; // UPDATE FOR KF
      u[1][0] = 0;
      u[2][0] = 0;
      u[3][0] = 0;
      
      powerSet(0, 0, 0, 0);
    } else {
      // Otherwise, motor power commands should be
      // chosen by the controller

      // FIXME
      // KF Filter update
      predict();
      update();
      
      o_x = x[0][0];
      o_y = x[1][0];
      o_z = x[2][0];
      psi = x[3][0];
      theta = x[4][0];
      phi = x[5][0];
      v_x = x[6][0];
      v_y = x[7][0];
      v_z =x[8][0];

      tau_x = 0.00303466f * (o_y - o_y_des) -0.00307415f * phi + 0.00141493f * v_y -0.00037938f * w_x;
      tau_y = -0.00341735f * (o_x - o_x_des) -0.00641297f * theta -0.00214120f * v_x -0.00100975f * w_y;
      tau_z = -0.00052902f * psi -0.00022111f * w_z;
      f_z = -0.47852401f * (o_z - o_z_des) -0.17814013f * v_z + 0.31019220f;

      tau_x = 0.00063277f * (o_y - o_y_des) -0.00368826f * phi + 0.00093606f * v_y -0.00073120f * w_x;
      tau_y = -0.00068347f * (o_x - o_x_des) -0.00398028f * theta -0.00101082f * v_x -0.00078839f * w_y;
      tau_z = -0.00013659f * psi -0.00016561f * w_z;
      f_z = -0.12789092f * (o_z - o_z_des) -0.15634548f * v_z + 0.31019220f;

      tau_x = 0.00002001f * (o_y - o_y_des) -0.00026050f * phi + 0.00003825f * v_y -0.00009941f * w_x;
      tau_y = -0.00002161f * (o_x - o_x_des) -0.00027992f * theta -0.00004124f * v_x -0.00010643f * w_y;
      tau_z = -0.00000432f * psi -0.00001720f * w_z;
      f_z = -0.00404427f * (o_z - o_z_des) -0.01649592f * v_z + 0.31019220f;

      tau_x = 0.00000200f * (o_y - o_y_des) -0.00007070f * phi + 0.00000573f * v_y -0.00005077f * w_x;
      tau_y = -0.00000216f * (o_x - o_x_des) -0.00007590f * theta -0.00000617f * v_x -0.00005431f * w_y;
      tau_z = -0.00000043f * psi -0.00000528f * w_z;
      f_z = -0.00040443f * (o_z - o_z_des) -0.00507341f * v_z + 0.31019220f;

      u[0][0] = tau_x; // UPDATE FOR KF
      u[1][0] = tau_y;
      u[2][0] = tau_z;
      u[3][0] = f_z;

        // FIXME
      m_1 = limitUint16( -3937708.6f * tau_x -3937708.6f * tau_y -48262548.3f * tau_z + 126262.6f * f_z );
      m_2 = limitUint16( -3937708.6f * tau_x + 3937708.6f * tau_y + 48262548.3f * tau_z + 126262.6f * f_z );
      m_3 = limitUint16( 3937708.6f * tau_x + 3937708.6f * tau_y -48262548.3f * tau_z + 126262.6f * f_z );
      m_4 = limitUint16( 3937708.6f * tau_x -3937708.6f * tau_y + 48262548.3f * tau_z + 126262.6f * f_z );
      // Apply motor power commands
      powerSet(m_1, m_2, m_3, m_4);