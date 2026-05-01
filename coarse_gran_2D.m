clear
clc

for flagg=1:1
    
Nfile = 371; % File averaging number
startfile = 300000; % Start file timestep
fileinc = 10000; % File timestep increment
dl = 0.0008; % Size of larger grains
ds = 0.0008; % Size of smaller grains
err = 0.1; % Grain size fluctuation
dllower = dl*(1-err); % Lower bound of larger grain size
dlupper = dl*(1+err); % Upper bound of larger grain size
dslower = ds*(1-err); % Lower bound of smaller grain size
dsupper = ds*(1+err); % Upper bound of smaller grain size
init_cl = 0.75; % Initial cl
init_cs = 1 - init_cl;
d = init_cl*dl + init_cs*ds; % Averaged grain size
exclude_top = 0; 
exclude_bottom = 0;
binsize = 4.0*dl; % Thickness of each bin, relative to the larger grain size
Nslice = 600; % Slice distance is approximately 0.6d now
area_bin = binsize*60*d;
%Vxl = zeros(Nfile,Nslice); % Final Vx of each slice
%Vxs = zeros(Nfile,Nslice);

Vx = zeros(Nfile,Nslice);
Vx_grad = zeros(Nfile,Nslice);
Vx_grad_grad = zeros(Nfile,Nslice);

dfield = zeros(Nfile,Nslice);

% --- Add this in the variables section ---
Sxx_all = zeros(Nfile, 1);
Syy_all = zeros(Nfile, 1);
Sxy_all = zeros(Nfile, 1);
time_all = zeros(Nfile, 1);
dt_val = 2e-07; % Timestep from your LAMMPS script

phi = zeros(Nfile,Nslice);
cl  = zeros(Nfile,Nslice);
cl_grad = zeros(Nfile,Nslice);

pressure = zeros(Nfile,Nslice);
pressure_grad = zeros(Nfile,Nslice);

z_slice = zeros(Nfile,Nslice);
z_spatial = zeros(Nfile,Nslice);

maxx=[];
   
% OUTSIDE the loop — keep only these declarations:
Sxx_bulk = zeros(Nfile, 1);
Syy_bulk = zeros(Nfile, 1);
Sxy_bulk = zeros(Nfile, 1);
time_vec = zeros(Nfile, 1);
dt_val = 2e-07;


for filecount = 1:Nfile
 
  currentfile = startfile + (filecount-1)*fileinc;
    textFileName = ['output-shear-' num2str(currentfile) '.txt'];
	if exist(textFileName, 'file')
        data = dlmread(textFileName,' ',9,0);
      
        % --- Inside the loop, after data = dlmread(...) ---

        % 1. Calculate the current time for this file
        time_vec(filecount) = currentfile * dt_val;

        % 2. Get current box dimensions
        Lx = 60 * d; % Ensure this matches your actual box width
        Ly = max(data(:,4)) - min(data(:,4));
        Area = Lx * Ly;





        % 1. Calculate the current time for this file
        time_vec(filecount) = currentfile * dt_val;

        % 2. Get current box dimensions
        Lx = 60 * d; 
        Ly = max(data(:,4)) - min(data(:,4));
        Area = Lx * Ly;

        % 3. Calculate Bulk Stresses (Sum / Area)
        % Column 11=Sxx, 12=Syy, 13=Sxy (c_peratom[4])
        Sxx_bulk(filecount) = sum(data(:,11)) / Area;
        Syy_bulk(filecount) = sum(data(:,12)) / Area;
        Sxy_bulk(filecount) = sum(data(:,13)) / Area;
        
        bottom = min(data(:,4)) + exclude_bottom;
        top = max(data(:,4)) - exclude_top;
        Delta = (top-bottom)/(Nslice-1); % Distance between two slices.
        y_slice_inst =top:-Delta:bottom; % Coords of all slices.
        
        Natoms = size(data,1); % Total grain number
        

    
    c=W/2;dummy = -c:del:c;
    lucy_fun = 1*(-3*abs(dummy/c).^4 + 8*abs(dummy/c).^3 -6*abs(dummy/c).^2 +1);
    norm = sum(abs(lucy_fun));
    lucy_fun = lucy_fun/norm;
    lucy_fun_grad = 1/norm*(+12/c*abs(dummy(1:n)/c).^3 - 24/c*abs(dummy(1:n)/c).^2 + 12/c*abs(dummy(1:n)/c).^1 );
    lucy_fun_grad = -[lucy_fun_grad, 1/norm*(-12/c*abs(dummy(n+1:2*n+1)/c).^3 + 24/c*abs(dummy(n+1:2*n+1)/c).^2 - 12/c*abs(dummy(n+1:2*n+1)/c).^1 )];
    
    lucy_fun_grad_grad = 1/norm*(+36/c^2 *abs(dummy(1:n)/c).^2 - 48/c^2*abs((dummy(1:n)/c)) + 12/c^2);
    
    lucy_fun_grad_grad = -[lucy_fun_grad_grad, -1/norm*(-36/c^2 *abs(dummy(n+1:2*n+1)/c).^2 + 48/c^2*abs((dummy(n+1:2*n+1)/c)) - 12/c^2)];
    
%     %%%%%%%%%%%Gaussian Function%%%%%%%%%%%%%%%%%%%%
    sigma=W/10; dummy = -W/2:del:W/2;
    gauss_fun = exp(-dummy.^2/2/sigma^2);
    norm = sum(abs(gauss_fun));
    gauss_fun = 1/norm*(gauss_fun);
    gauss_fun_grad = 1/norm/sigma^2* exp(-dummy.^2/2/sigma^2).*dummy;
    
    gauss_fun_grad_grad = 1/norm/sigma^4 *exp(-dummy.^2/2/sigma^2).*(dummy.^2 - sigma^2);

        phi_inst = zeros(Nslice,1);
        cl_inst =zeros(Nslice,1);
        cl_grad_inst = zeros(Nslice,1);
      	Vx_inst = zeros(Nslice,1);	
		Vx_grad_inst = zeros(Nslice,1);
        Vx_grad_grad_inst = zeros(Nslice,1);
        dfield_inst = zeros(Nslice,1);
      
        pressure_inst = zeros(Nslice,1);
        pressure_grad_inst = zeros(Nslice,1);
   for k=1:Nslice
       
       cross_sec_inst=zeros(m,1);
       cross_sec_l = zeros(m,1);
	   cross_sec_s = zeros(m,1);	
       
       Vx_sub_inst = zeros(m,1);
       pressure_sub_inst = zeros(m,1);
	
	
       Zsub =y_slice_inst(k)-W/2 :del: y_slice_inst(k)+W/2;
     for m=1:2*n+1
       for line=1:Natoms
           
           
           Zcoord = data(line,4);
           diam = data(line,9);
           radius = diam/2;
           grain_area = radius*radius*pi;
           h=abs(Zcoord-Zsub(m));
            

               if h < (radius + binsize/2) % This grain intercept a bin
                    if h > binsize/2 % Center of grain outside the bin
                        arc_height = (radius + binsize/2 - h);
                        cross_sec = radius*radius*acos((radius-arc_height)/radius) - (radius - arc_height)*sqrt(2*radius*arc_height - arc_height*arc_height);
                    elseif h > (binsize/2 - radius) % Center of grain inside the bin but one segment is outside
                        arc_height = (radius - binsize/2 + h);
                        cross_sec = radius*radius*pi - radius*radius*acos((radius-arc_height)/radius) + (radius - arc_height)*sqrt(2*radius*arc_height - arc_height*arc_height);
                    else % Fully inside the bin
                        cross_sec = grain_area;
                    end
                    cross_sec_inst(m) = cross_sec_inst(m) + cross_sec; % Add on to total cross section of this slice
                    Vx_sub_inst(m) = Vx_sub_inst(m) + cross_sec*data(line,6);
                    pressure_sub_inst(m) = pressure_sub_inst(m) - cross_sec*(data(line,11)+data(line,12))/2/grain_area;
                   
                    if (diam > dllower) && (diam < dlupper)
                        cross_sec_l(m) = cross_sec_l(m) + cross_sec;
                       % Vxl_inst(slice) = Vxl_inst(slice) + cross_sec*data(line,6);
                        % This elseif is actually unneccesary
                    elseif (diam > dslower) && (diam < dsupper)
                        cross_sec_s(m) = cross_sec_s(m) + cross_sec;
                       % Vxs_inst(slice) = Vxs_inst(slice) + cross_sec*data(line,6);
                    else
                        fprintf('Wrong grain size.\n');
                    end
               end
           
       end
     end
     
    if (abs(flagg-1)<1e-3)   
    %%%%%%%%%%%Lucy fun%%%%%%%%%%%%%     
    wei_fun = lucy_fun;
    wei_fun_grad = -lucy_fun_grad;
    wei_fun_grad_grad = -lucy_fun_grad_grad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
 elseif(abs(flagg-2)<1e-3)    
%%%%%%%%%%%Gauss fun%%%%%%%%%%%%%     
   wei_fun = gauss_fun;
   wei_fun_grad = -gauss_fun_grad;
   wei_fun_grad_grad = -gauss_fun_grad_grad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end

         %%  vol_frac_inst(k) =wei_fun * (cross_sec_inst./area);
           phi_inst(k) = wei_fun * (cross_sec_inst./area_bin);
           cl_inst(k) = wei_fun * (cross_sec_l./cross_sec_inst);
           
           cl_grad_inst(k) = wei_fun_grad * (cross_sec_l./cross_sec_inst);
           
		  % Vxs_inst(k) = wei_fun * ( Vxs_sub_inst./cross_sec_s);

          % Vxl_inst(k) = wei_fun * ( Vxl_sub_inst./cross_sec_l);		   
		   
		    Vx_inst(k) = wei_fun * (Vx_sub_inst./cross_sec_inst);
           
            Vx_grad_inst(k) = wei_fun_grad * (Vx_sub_inst./cross_sec_inst);
            
            Vx_grad_grad_inst(k) = wei_fun_grad_grad * (Vx_sub_inst./cross_sec_inst);
            
            dfield_inst(k) = (1-cl_inst(k))*ds + cl_inst(k)*dl;
            
            pressure_inst(k) = wei_fun * (pressure_sub_inst./area_bin);
            
            pressure_grad_inst(k) = wei_fun_grad * (pressure_sub_inst./area_bin);
     
   end
   
   filecount

        Vx(filecount,:) = Vx_inst';
        Vx_grad(filecount,:) = Vx_grad_inst';
        Vx_grad_grad(filecount,:) = Vx_grad_grad_inst';

        phi(filecount,:) = phi_inst';
        cl(filecount,:) = cl_inst';
        cl_grad(filecount, :) = cl_grad_inst';
        
        dfield(filecount,:) = dfield_inst';
        
        pressure(filecount,:) = pressure_inst';
        pressure_grad(filecount,:) = pressure_grad_inst';
        
        z_slice(filecount,:) = y_slice_inst';
        
        z_spatial(filecount,:) = y_slice_inst' - y_slice_inst(end) + exclude_bottom;
    else
		fprintf('File %s does not exist.\n', textFileName);
    
    end 

end


phi_avg = mean(phi,1);
cl_avg = mean(cl,1);
cl_grad_avg = mean(cl_grad,1);
dfield_avg = mean(dfield,1);
z_avgg       = mean(z_slice, 1);
z_avg        = z_avgg(1) - z_avgg;
z_norm       = mean(z_spatial, 1) ./ d;


% ==========================================================
% FORCE CONSTANT GLOBAL GAMMA DOT (SLOPE BETWEEN ENDPOINTS)
% ================================================================
% CONSTANT GAMMA DOT via linear regression on physical height (y)
% ================================================================
% Use 10% trim to avoid wall effects
trim = round(0.1 * Nslice); 

% Use physical distance in meters for the fit to get 1/s units
z_phys = mean(z_spatial, 1); 
Vx_avg_temp = mean(Vx, 1);

% Extract the middle portion for regression
z_fit  = z_phys(trim:end-trim);
Vx_fit = Vx_avg_temp(trim:end-trim);

% Perform linear fit: Vx = p(1)*y + p(2)
% p(1) is dVx/dy, which is exactly gamma_dot
p = polyfit(z_fit, Vx_fit, 1);
gamma_dot_global = p(1);

% Display the value in the command window to verify it isn't tiny
fprintf('Calculated Global Shear Rate: %.4f s^-1\n', gamma_dot_global);

% Fill the averaging arrays with this constant value
Vx_grad_avg = ones(1, Nslice) * gamma_dot_global;

for f = 1:Nfile
    Vx_grad(f, :) = gamma_dot_global;
end

Vx_grad_grad_avg = zeros(1, Nslice);
Vx_grad_grad     = zeros(Nfile, Nslice);

% ==========================================================

% 3. Final averages that DO depend on Vx_grad
Vx_avg           = mean(Vx, 1);
Vx_grad_avg      = mean(Vx_grad, 1);
Vx_grad_grad_avg = mean(Vx_grad_grad, 1);
pressure_avg     = mean(pressure, 1);
pressure_grad_avg = mean(pressure_grad, 1);



if (abs(flagg-1)<1e-3)
%%%%%%  writing fields using lucy kernel %%%%%%%%%%%%%%
dlmwrite('dia_lucy_transient.txt',dfield);
dlmwrite('cl_lucy_transient.txt',cl);
dlmwrite('cl_grad_lucy_transient.txt',cl_grad);
dlmwrite('Vx_lucy_transient.txt', Vx);
dlmwrite('Vx_grad_lucy_transient.txt',Vx_grad);
dlmwrite('Vx_grad_grad_lucy_transient.txt',Vx_grad_grad);
dlmwrite('pressure_lucy_transient.txt',pressure);
dlmwrite('pressure_grad_lucy_transient.txt',pressure_grad);
dlmwrite('z_norm_transient.txt',z_norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(abs(flagg-2)<1e-3)
%%%%%%  writing fields using gauss kernel %%%%%%%%%%%%%%
dlmwrite('dia_gauss.txt',dfield);
dlmwrite('cl_gauss.txt',cl_avg);
dlmwrite('cl_grad_gauss.txt',cl_grad);
dlmwrite('Vx_gauss.txt',Vx);
dlmwrite('Vx_grad_gauss.txt',Vx_grad);
dlmwrite('Vx_grad_grad_gauss.txt',Vx_grad_grad);
dlmwrite('pressure_gauss.txt',pressure);
dlmwrite('pressure_grad_gauss.txt',pressure_grad);
dlmwrite('z_norm.txt',z_norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end


%% VISUALIZATION SECTION

% 1. KINEMATICS AND SPECIES DISTRIBUTION
figure('Name', 'Kinematics and Concentration', 'NumberTitle', 'off');

% Velocity Profile
subplot(2,2,1)
plot(Vx_avg, z_norm, 'b-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$V_x$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Velocity Profile')

% Shear Rate (Velocity Gradient)
subplot(2,2,2)
plot(Vx_grad_avg, z_norm, 'b-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$\dot{\gamma}$ (1/s)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Local Shear Rate')

% Large Grain Concentration
subplot(2,2,3)
plot(cl_avg, z_norm, 'g-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$c^l$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Concentration $c^l$')

% Concentration Gradient
subplot(2,2,4)
plot(cl_grad_avg, z_norm, 'g-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$\nabla c^l$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Concentration Gradient')

% 2. FORCES AND GEOMETRY
figure('Name', 'Forces and Scale', 'NumberTitle', 'off');

% Velocity Second Gradient (Curvature)
subplot(2,2,1)
plot(Vx_grad_grad_avg, z_norm, 'r-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$\nabla^2 V_x$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Velocity Curvature')

% Mean Diameter Field
subplot(2,2,2)
plot(dfield_avg, z_norm, 'r-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$d_{field}$ (m)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Mean Diameter')

% Pressure Profile
subplot(2,2,3)
plot(pressure_avg, z_norm, 'k-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('Pressure (Pa)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Normal Stress (Pressure)')

% Pressure Gradient
subplot(2,2,4)
plot(pressure_grad_avg, z_norm, 'k-', 'LineWidth', 2)
grid on; set(gca,'FontSize', 12)
xlabel('$\nabla P$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Height ($y/d$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Pressure Gradient')



% Stress Visualization Section
figure('Name', 'Stress Evolution vs Time', 'NumberTitle', 'off');

subplot(3,1,1)
plot(time_vec, Sxx_bulk, 'r', 'LineWidth', 1.5); hold on;
plot(time_vec, Syy_bulk, 'b', 'LineWidth', 1.5);
grid on;
ylabel('Normal Stress (Pa)');
legend('\sigma_{xx}', '\sigma_{yy}');
title('Normal Stress vs Time');

subplot(3,1,2)
plot(time_vec, Sxy_bulk, 'm', 'LineWidth', 1.5);
grid on;
ylabel('Shear Stress \sigma_{xy} (Pa)');
title('Shear Stress vs Time');

subplot(3,1,3)
% Effective Friction Coefficient (Shear / Normal)
plot(time_vec, Sxy_bulk ./ Syy_bulk, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('\mu = \sigma_{xy} / \sigma_{yy}');
title('Friction Coefficient (\mu) Evolution');



% Export Stress Data to Text File
% Columns: Time | Sxx | Syy | Sxy
output_matrix = [time_vec, Sxx_bulk, Syy_bulk, Sxy_bulk];

% Save as a tab-delimited file
dlmwrite('stress_history_output.txt', output_matrix, 'delimiter', '\t', 'precision', '%.10f');

fprintf('Success! Stress data saved to: stress_history_output.txt\n');