
 function  ss = SEIR_ss(params,data)

% N = 1e4; E0 = 0; I0 = 1; R0 = 0; S0 = N - E0 - I0 - R0; 

  N =  data.N; 
  Y0 = data.Y0;
  range = data.range;

  t_data = data.xdata;
  ydata = data.ydata;

  ode_options = odeset('RelTol',1e-6);
  [t,Y] = ode45(@(t,y) SEIR_rhs(t,y,params,N),range,Y0,ode_options); 
    s = spline(t,Y(:,3),t_data);
  ss = sum((s - ydata).^2);