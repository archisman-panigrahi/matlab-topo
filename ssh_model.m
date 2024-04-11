clear all;
global t1 = 1;
global t = 1;

function energy_with_m(nsites, n_points_in_plot)
global t;
global t1;
no_of_points = n_points_in_plot;
m_in_paper = linspace(-2,2,no_of_points);
m = nsites; n = 1; p = 1;
lattices1D;
angmom;
for i=1:no_of_points
	h_ssh = kron(m_in_paper(i) * M1D + t1 * CX1Dp, sigma_x) + t * kron(SX1Dp, sigma_y);
	energy_eigenvalues = eig(h_ssh);

	for r=1:2*m
		energies( (i-1)*(2*m) + r ) = energy_eigenvalues(r);
	end

end


xvar = kron(m_in_paper,linspace(1,1,2*m));
graphics_toolkit('gnuplot');
text("interpreter","latex");
figure(1);

scatter(xvar,energies,'.');


xlabel('$m$', "fontsize", 20);
ylabel('$E/t_1$', "fontsize", 20);
box on;
set(gca, "linewidth", 2, "fontsize", 20);
end


function energy_eigenvalues_in_momentum_space(nsites, n_points_in_plot)
	global t;
	global t1;
	no_of_points = n_points_in_plot;
	m_in_paper = linspace(-2,2,no_of_points);
	k_x = linspace(-pi,pi,nsites);
	% h_ssh = kron(m_in_paper(i) * M1D + t1 * CX1Dnp, sigma_x) + t * kron(SX1Dnp, sigma_y);
	for l = 1 : no_of_points
		for k = 1 : nsites
			en_eig = hypot(m_in_paper(l) + t1 * cos(k_x(k)), t * sin(k_x(k)));
			eigenvalues((l-1)*(2*nsites) + (2*k - 1)) = en_eig;
			eigenvalues((l-1)*(2*nsites) + 2*k) = -en_eig;
		end
	end
	xvar = kron(m_in_paper, linspace(1,1,2*nsites));
	figure;
	scatter(xvar, eigenvalues,'.');

	text("interpreter","latex");
	xlabel('$m$');
	ylabel('allowed energies in units of $t_1$');
end

function print_energy_with_m(nsites, n_points_in_plot, name)
	% graphics_toolkit('gnuplot'); %Always use gnuplot for printing
	energy_with_m(nsites, n_points_in_plot);
	cd saved_plots/ssh_model;
	print(figure(1),'-dpdflatexstandalone', name);
	system (strcat("pdflatex\t", name));
	cd ../../;
end


function Energyspectra2D(m_in_paper, nsites, n_points_in_plot, t, t1, t2, name)
%m_in_paper is specified
no_of_points = n_points_in_plot;
m = nsites; n = 1; p = 1;
lattices;
angmom;
ky = linspace(-pi ,pi ,no_of_points);
for i=1:no_of_points
	h_ssh = kron(m_in_paper * M1D + t1 * CX1Dnp + t2 * cos(ky(i)) * M1D, sigma_x) + t * kron(SX1Dnp, sigma_y);
	energy_eigenvalues = eig(h_ssh);
	for r=1:2*m
		energies( (i-1)*(2*m) + r ) = energy_eigenvalues(r);
	end

end


xvar = kron(ky/pi, linspace(1,1,2*m));

figure(1);

scatter(xvar,energies,'.');

text("interpreter","latex");
xlabel('$k_y(\pi)$');
ylabel('allowed energies in units of $t_1$');
cd saved_plots/ssh_model;

	print(figure(1),'-dpdflatexstandalone', name);
	system (strcat("pdflatex\t", name));
	cd ../../
end


function ImportantBands2D(nsites,m_in_paper,namesave)
close all;
t = 1;
t1 = 1;
t2 = 0.5;
% no_of_points = n_points_in_plot;
m = nsites; n = 1; p = 1;
lattices;
angmom;
% m_in_paper = [-1.75,-1,0,1,1.75];
ky = linspace(-pi ,pi , 100);
for i=1:100
	h_ssh = kron(m_in_paper * M1D + t1 * CX1Dnp + t2 * cos(ky(i)) * M1D, sigma_x) + t * kron(SX1Dnp, sigma_y);
	energy_eigenvalues = eig(h_ssh);
	energy_to_plot1(i) = energy_eigenvalues(m-1);
	energy_to_plot2(i) = energy_eigenvalues(m);
	energy_to_plot3(i) = energy_eigenvalues(m+1);
	energy_to_plot4(i) = energy_eigenvalues(m+2);	
end



figure(1);

plot(ky/pi,energy_to_plot1,'*');
hold on
plot(ky/pi,energy_to_plot2,'+');
hold on
plot(ky/pi,energy_to_plot3,'--','color','k');
hold on
plot(ky/pi,energy_to_plot4,'*');
text("interpreter","latex");
xlabel('$k_y(\pi)$');
ylabel('allowed energies in units of $t_1$');
legend('1','2','3','4');
title(strcat('$m = ',num2str(m_in_paper),'$'));

text("interpreter","latex");

cd saved_plots/ssh_model/spectral_gap/bands;
print(figure(1),'-dpdflatexstandalone', namesave);
system (strcat("pdflatex\t", namesave));

cd ../../../..
end

function spectralGap2D(nsites)
close all;
t = 1;
t1 = 1;
t2 = 0.5;
% no_of_points = n_points_in_plot;
m = nsites; n = 1; p = 1;
lattices;
angmom;
m_in_paper = [-1.75,-1,-0.5,0,0.5,1,1.75];
ky = linspace(-pi , pi, 100);
namesave = ['spectgapm-1pt75';'spectgapm-1';'spectgapm-0pt5';'spectgapm0';'spectgapm0pt5';'spectgapm1';'spectgapm1pt75'];
for k = 1:7
for i=1:100
	h_ssh = kron(m_in_paper(k) * M1D + t1 * CX1Dnp + t2 * cos(ky(i)) * M1D, sigma_x) + t * kron(SX1Dnp, sigma_y);
	energy_eigenvalues = eig(h_ssh);
	energy_diff(i) = energy_eigenvalues(m+2) - energy_eigenvalues(m-1);
	%second excited state(first excited non-edge state) - state before edge state (last occupied non edge state)
end



figure(k);

plot(ky/pi,energy_diff);
text("interpreter","latex");
xlabel('$k_y(\pi)$');
ylabel('allowed energies in units of $t_1$');
title(strcat('$m = ',num2str(m_in_paper(k)),'$'));
grid on;

text("interpreter","latex");
% dirname = strcat('Lx',num2str(m));
dirname = strcat('Lx',num2str(m));
cd saved_plots/ssh_model/spectral_gap;
mkdir(dirname);
cd(dirname);
% Change the name in the above two lines accordingly and call ssh_model once again before plotting
print(figure(k),'-dpdflatexstandalone', namesave(k,:));
system (strcat("pdflatex\t", namesave(k,:)));

cd ../../../..
end
end

function Energyspectrawithm2D(no_of_points, nsitesx, nsitesy, t, t1, t2)
%Solved in 2D real space for different values of m_in_paper
%Diagonalizing the whole 2d eigenspace to obtain eigenstates
% no_of_points = n_points_in_plot;
	m = nsitesx; n = nsitesy; p = 1;
	lattices;
	angmom;
	m_in_paper = linspace(-2,2,no_of_points);
	
	for i=1:no_of_points
		h_ssh = kron(m_in_paper(i) * M2D + t1 * CX2Dnp + t2 * CY2Dnp, sigma_x) + t * kron(SX2Dnp, sigma_y);
		%There are zero energy states for both periodic and non-periodic. Check if they are edge states.
		energy_eigenvalues = eig(h_ssh);
		for r=1:2*m*n
			energies( (i-1)*(2*m*n) + r ) = energy_eigenvalues(r);
		end
	end
	% figure(1);
	% scatter(linspace(1,2*m*n,2*m*n),diag(energy_eigenvalues));


	xvar = kron(m_in_paper,linspace(1,1,2*m*n));

	figure;
	scatter(xvar,energies,'.');

	text("interpreter","latex");
	xlabel('$m$');
	ylabel('allowed energies in units of $t_1$');
	title('Energy Spectra as a function of $m$ for 2D open lattice')

% 	cd saved_plots;

% 	print(figure(1),'-dpdflatexstandalone', name);
% 	system (strcat("pdflatex\t", name));
% 	cd ..
end

function showeigenstates2D(m_in_paper, nsitesx, nsitesy, t, t1, t2, name1, name2, name3, name4)
%Diagonalizing the whole 2d eigenspace to obtain eigenstates
% no_of_points = n_points_in_plot;
	close all;
	m = nsitesx; n = nsitesy; p = 1;
	lattices;
	angmom;

	h_ssh = kron(m_in_paper * M2D + t1 * CX2Dnp + t2 * CY2Dnp, sigma_x) + t * kron(SX2Dnp, sigma_y);
	[eigenstates, energy_eigenvalues] = eig(h_ssh);

	state_of_interest = eigenstates(:, m*n);

	for k=1:m*n
		prob_density_left(k) = state_of_interest(2*k - 1) * conj(state_of_interest(2*k - 1));
		prob_density_right(k) = state_of_interest(2*k) * conj(state_of_interest(2*k));
	end

	for k=1:m
		for l=1:n
			prob_density_left_z(k,l) = prob_density_left((l-1)*m + k);
			prob_density_right_z(k,l) = prob_density_right((l-1)*m + k);
		end
	end
	[xx, yy] = meshgrid(linspace(1,m,m),linspace(1,n,n));
	figure;
	surf(xx,yy,prob_density_left_z);
	xlabel('$x$');
	ylabel('$y$');
	zlabel('probability density');
	title('eigenstate of orbital A')
	colorbar ("peer", gca, "northoutside");
	% legend(strcat('$m = $',num2str(m_in_paper)));

	figure;
	surf(xx,yy,prob_density_right_z);
	xlabel('$x$');
	ylabel('$y$');
	zlabel('probability density');
	title('eigenstate of orbital B');
	colorbar ("peer", gca, "northoutside");
	% legend(strcat('$m = $',num2str(m_in_paper)));

	figure;
	contourf(xx,yy,prob_density_left_z);
	xlabel('$x$');
	ylabel('$y$');
	zlabel('probability density');
	title('eigenstate of orbital A')
	colorbar ("peer", gca, "northoutside");
	axis equal;
	% legend(strcat('$m = $',num2str(m_in_paper)));

	figure;
	contourf(xx,yy,prob_density_right_z);
	xlabel('$x$');
	ylabel('$y$');
	zlabel('probability density');
	title('eigenstate of orbital B');
	colorbar ("peer", gca, "northoutside");
	axis equal;
	% legend(strcat('$m = $',num2str(m_in_paper)));
	
	cd saved_plots/eigenstates2d;

	print(figure(1),'-dpdflatexstandalone', name1);
	system (strcat("pdflatex\t", name1));
	print(figure(2),'-dpdflatexstandalone', name2);
	system (strcat("pdflatex\t", name2));
	print(figure(3),'-dpdflatexstandalone', name3);
	system (strcat("pdflatex\t", name3));
	print(figure(4),'-dpdflatexstandalone', name4);
	system (strcat("pdflatex\t", name4));
	cd ../..
end

% Percentage of Edge States in 2D

function edgestates2D(nsitesx, nsitesy, m_in_paper)
	t1 = 1;
	t = 1;
	t2 = 0.5;
	m = nsitesx;
	n = nsitesy;
	p = 1;
	lattices;
	angmom;
	close all;
	h_ssh_real2D = kron(m_in_paper * M2D + t1 * CX2Dnp + t2 * CY2Dnp, sigma_x) + t * kron(SX2Dnp, sigma_y);
	[eigenstates, energy_eigenvalues] = eig(h_ssh_real2D);

	state_of_interest = eigenstates(:, m*n);
	for k=1:n
		perc_eigfunc_loc_at_left_edgeA(k) = abs(state_of_interest(1 + 2*m*(k-1)))^2 * 100;
		perc_eigfunc_loc_at_right_edgeA(k) = abs(state_of_interest(2*m*k - 1))^2 * 100;
		perc_eigfunc_loc_at_left_edgeB(k) = abs(state_of_interest(2 + 2*m*(k-1)))^2 * 100;
		perc_eigfunc_loc_at_right_edgeB(k) = abs(state_of_interest(2*m*k))^2 * 100;
	end
	graphics_toolkit('fltk');
	figure;
	xvar = linspace(1,n,n);
	plot(xvar,perc_eigfunc_loc_at_left_edgeA,'--','color','blue', xvar, perc_eigfunc_loc_at_left_edgeB, 'x','color','red', xvar, perc_eigfunc_loc_at_right_edgeA, '--','color','green', xvar, perc_eigfunc_loc_at_right_edgeB,'*','color','black');
	title(strcat('percentage of occupied state with highest energy localized at edges for $m = ', num2str(m_in_paper), ', L_x =', num2str(nsitesx), ' $'))
	ylabel('percentage of occupied state with highest energy localized at edges');
	xlabel('$y$');
	legend('leftA','leftB','rightA','rightB');
end

%Variation of k_y near zero - eigenstate at leftmost site for different m_in_paper
function variationOfPsiWithKy(nsitesx)
	t1 = 1;
	t = 1;
	t2 = 0.5;
	m = nsitesx;
	n = 1;
	p = 1;
	lattices;
	angmom;
	close all;
	m_in_paper = [-1.75, -1, 0, 1, 1.75];
	k_y = linspace(-pi,pi,20);
	for k=1:20
		for i= 1:5
			h_ssh = kron(m_in_paper(i) * M1D + t1 * CX1Dnp + t2 * cos(k_y(k)), sigma_x) + t * kron(SX2Dnp, sigma_y);
			[eigenstates, energy_eigenvalues] = eig(h_ssh);
			to_plotA(k,i) = abs(eigenstates(1,m))^2;
			to_plotB(k,i) = abs(eigenstates(2,m))^2;
		end
	end
	graphics_toolkit('gnuplot')
	figure(1);
	plot(k_y/pi, to_plotA(:,1),'color','red',k_y/pi, to_plotA(:,2),'color','blue',k_y/pi, to_plotA(:,3),'color','black',k_y/pi, to_plotA(:,4),'color','green',k_y/pi, to_plotA(:,5),'color','cyan')
	legend('-1.75','-1','0','1','1.75');
	xlabel('$k_y(\pi)$');
	ylabel('probability density of orbital A for 1st site')
	title(strcat('$L_x = ',num2str(m),'$'))	
	figure(2);
	plot(k_y/pi, to_plotB(:,1),'color','red',k_y/pi, to_plotB(:,2),'color','blue',k_y/pi, to_plotB(:,3),'color','black',k_y/pi, to_plotB(:,4),'color','green',k_y/pi, to_plotB(:,5),'color','cyan')
	legend('-1.75','-1','0','1','1.75');
	xlabel('$k_y(\pi)$');
	ylabel('probability density of orbital B for 1st site')	
	title(strcat('$L_x = ',num2str(m),'$'))
	cd saved_plots/ssh_model/statevsky;
	name1 = strcat('siteAL_x',num2str(m))
	name2 = strcat('siteBL_x',num2str(m))
	print(figure(1),'-dpdflatexstandalone',name1);
	print(figure(2),'-dpdflatexstandalone',name2);
	system(strcat("pdflatex\t",name1))
	system(strcat("pdflatex\t",name2))
	cd ../../..
end


%Visualize eigenstates

function showeigenstates_np(no_of_sites, m_in_paper,t1,t,name)
close all;
% global t;
% global t1;
% t1 = 1;
% t = 2;
m = no_of_sites;
lattices;
angmom;


	h_ssh = kron(m_in_paper * M1D + t1 * CX1Dnp, sigma_x) + t * kron(SX1Dnp, sigma_y);
	[eigenstates,energy_eigenvalues] = eig(h_ssh);
	lowest_energy_eigstate = eigenstates(:, m); %octave lists them in increasing order
	first_excited_state = eigenstates(:, m + 1);
	second_lowest_energy_state = eigenstates(:, m - 1);

	for k=1:m
		%combining the two orbitals at every site
		prob_density_left(k) = lowest_energy_eigstate(2*k - 1) * conj(lowest_energy_eigstate(2*k - 1));
		prob_density_right(k) = lowest_energy_eigstate(2*k) * conj(lowest_energy_eigstate(2*k));
	end
	xvar = linspace(1,m,m);
	prob_density_total = prob_density_right + prob_density_left;
	graphics_toolkit('gnuplot')

	figure;
	plot(xvar,prob_density_left,'-','color','blue',xvar,prob_density_right,'-','color','red', xvar, prob_density_total,'--', 'color','black');
	xlabel("site number",'fontsize',20);
	ylabel("probablity density",'fontsize',20);
	axis([1 m]);
	legend({'left','right','total'});
	title(strcat('lowest absolute energy eigenstate for non-periodic lattice, nsites = ', num2str(m), ' & m = ',num2str(m_in_paper', '%0.2f')),'fontsize',20);
	grid on;
	set(gca, "linewidth", 2, "fontsize", 20);

	for k=1:m
		%combining the two orbitals at every site
		prob_density_left(k) = first_excited_state(2*k - 1) * conj(first_excited_state(2*k - 1));
		prob_density_right(k) = first_excited_state(2*k) * conj(first_excited_state(2*k));
	end
	xvar = linspace(1,m,m);
	prob_density_total = prob_density_right + prob_density_left;
	figure;
	plot(xvar,prob_density_left,'-','color','blue',xvar,prob_density_right,'-','color','red', xvar, prob_density_total,'--', 'color','black');
	xlabel("site number",'fontsize',20);
	ylabel("probablity density",'fontsize',20);
	axis([1 m]);
	legend({'left','right','total'});
	title(strcat('first excited state for non-periodic lattice, nsites = ', num2str(m), ' & m = ',num2str(m_in_paper', '%0.2f')),'fontsize',20);
	grid on;
	set(gca, "linewidth", 2, "fontsize", 20);

	for k=1:m
		%combining the two orbitals at every site
		prob_density_left(k) = second_lowest_energy_state(2*k - 1) * conj(second_lowest_energy_state(2*k - 1));
		prob_density_right(k) = second_lowest_energy_state(2*k) * conj(second_lowest_energy_state(2*k));
	end
	xvar = linspace(1,m,m);
	prob_density_total = prob_density_right + prob_density_left;
	figure;
	plot(xvar,prob_density_left,'-','color','blue',xvar,prob_density_right,'-','color','red', xvar, prob_density_total,'--', 'color','black');
	xlabel("site number",'fontsize',20);
	ylabel("probablity density",'fontsize',20);
	axis([1 m]);
	legend({'left','right','total'});
	title(strcat('second lowest abs energy eigenstate for non-periodic lattice, nsites = ', num2str(m), ' & m = ',num2str(m_in_paper', '%0.2f')),'fontsize',20);
	set(gca, "linewidth", 2, "fontsize", 20);
	grid on;
cd saved_plots/ssh_model/eigenstates1d
name1 = strcat('lowestenergystatelx',num2str(m),name);
name2 = strcat('firstexstatelx',num2str(m),name);
name3 = strcat('seclowenergystatelx',num2str(m),name);
print(figure(1),'-dpdflatexstandalone',name1);
% print(figure(2),'-dpdflatexstandalone',name2);
% print(figure(3),'-dpdflatexstandalone',name3);
system(strcat("pdflatex\t",name1));
% system(strcat("pdflatex\t",name2));
% system(strcat("pdflatex\t",name3));
cd ../../..
end

function showeigenstates_p(no_of_sites, m_in_paper)
global t;
global t1;

m = no_of_sites;
lattices;
angmom;
figure;

	h_ssh = kron(m_in_paper * M1D + t1 * CX1Dp, sigma_x) + t * kron(SX1Dp, sigma_y);
	[eigenstates,energy_eigenvalues] = eig(h_ssh);
	energy_eigenvalues = diag(energy_eigenvalues);
	[min_energy, idx] = min(abs(energy_eigenvalues));
	lowest_energy_eigstate = eigenstates(:,idx);

	prob_density = lowest_energy_eigstate .* conj(lowest_energy_eigstate);
	for k=1:m
		%combining the two orbitals at every site
		prob_density_combined(k) = prob_density(2*k) + prob_density(2*k - 1);
	end
	plot(linspace(1,m,m),prob_density_combined,'--');
	xlabel("site number");
	ylabel("probablity density");
	axis([1 m]);

legend(strcat('m = ',num2str(m_in_paper', '%0.2f')));
title(strcat('probability density of lowest absolute energy eigenstate for periodic lattice with nsites = ', num2str(m)));
grid on;
end

%%Bandstructure
%2D
function bandstructure2D()
	close all;
	t = 1;
	t1 = 1;
	t2 = 0.5;
	% Gamma(0,0) --> X(pi,0) --> M(pi,pi) --> Y(0,pi) --> Gamma(0,0) --> M(pi,pi)
	% h_ssh = kron(m_in_paper * M1D + t1 * CX1Dnp + t2 * cos(ky(i)) * M1D, sigma_x) + t * kron(SX1Dnp, sigma_y);

	x1 = linspace(0,pi,50);
	y1 = linspace(0,0,50);

	x2 = linspace(pi,pi,50);
	y2 = linspace(0,pi,50);

	x3 = linspace(pi,0,50);
	y3 = linspace(pi,pi,50);

	x4 = linspace(0,0,50);
	y4 = linspace(pi,0,50);

	x5 = linspace(0,pi,50);
	y5 = linspace(0,pi,50);

	k_x = [x1,x2,x3,x4,x5];
	k_y = [y1,y2,y3,y4,y5];

	m_in_paper(1) = 1;
	m_in_paper(2) = -1;

	for p=1:length(k_x)
		E1 = hypot(m_in_paper(1) + t1 * cos(k_x) + t2 * cos(k_y), t * sin(k_x));
		E2 = - hypot(m_in_paper(1) + t1 * cos(k_x) + t2 * cos(k_y), t * sin(k_x));
		E3 = hypot(m_in_paper(2) + t1 * cos(k_x) + t2 * cos(k_y), t * sin(k_x));
		E4 = - hypot(m_in_paper(2) + t1 * cos(k_x) + t2 * cos(k_y), t * sin(k_x)); 
	end
	% graphics_toolkit('gnuplot');
	text("interpreter","latex");
	figure(1);
	plot(linspace(1,length(k_x),length(k_x)), E1,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E2, 'linewidth', 4, 'color', 'blue');
	hold on
	plot(linspace(1,length(k_x),length(k_x)), E3,'linewidth', 4, 'color', 'red', linspace(1,length(k_x),length(k_x)), E4, 'linewidth', 4, 'color', 'red');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150 200 250]);
	set(gca, 'xticklabel',({'','','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-3.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-3.35, '$X$', 'fontsize', 30)
	text(100,-3.35, '$M$', 'fontsize', 30)
	text(150,-3.35, '$Y$', 'fontsize', 30)
	text(200,-3.35, '$\Gamma$', 'fontsize', 30)
	text(245, -3.35, '$M$', 'fontsize', 30)

	cd saved_plots/ssh_model/bandstructure
	print(figure(1),"-dpdflatexstandalone","bands")
	system("pdflatex bands")
	cd ../../..
end

%%%%%%%%Polarization%%%%%%%
%%%2D%%%

%Note, n1,px1 etc is the method in paper. n, px are my own method which gives different results for periodic and open
%boundary conditions
function polarization2D(nsites,m_in_paper,val_m_in_paper)
close all;
t = 1;
t1 = 1;
t2 =0.5;
no_of_points = 30;

k_y = linspace(-pi,pi,no_of_points);
m = nsites; n = 1; p = 1;
lattices;
angmom;
for i=1:no_of_points
	h_ssh = kron(m_in_paper * M1D + t1 * CX1Dnp + t2 * cos(k_y(i)) * M1D, sigma_x) + t * kron(SX1Dnp, sigma_y);
	[eigenstates, energy_eigenvalues] = eig(h_ssh);
	xmatrix = kron(diag(linspace(1,m,m)), eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%My method which gives different results for periodic and open boundary conditions
	%define number_at eigenstates
	% for k=1:m
	% 	number_at_diagonal(k) = 1;
	% 	number_at_diagonal(k+m) = 0;
	% end
	% number_at_diag = diag(number_at_diagonal);
	% eigenstates

	% number_at = eigenstates * number_at_diag * ctranspose(eigenstates);


	% n(i) = trace(number_at * xmatrix)/m;

	%average value of number_at is (0 + 1)/2 = 1/2; Then n_al = trace(1/2 * xmatrix)/m = 1/2 * 2 * m * (m+1)/m = (m+1)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%New method, according to BR's papers
	filled_eigenstates = eigenstates(:,1:m);
	x_exp = transpose(filled_eigenstates) * expm((2 * pi * j/m) * xmatrix) * conj(filled_eigenstates);
	n1(i) = real((-j/(2*pi)) * trace(logm(diag(eig(x_exp)))));

end
n_al = (m+1)/2;
% n = round(n * 100)/100; %This will round it of to 2 decimal places, otherwise values like -2e-5, 
n1 = round(n1 * 100)/100; %which should really be 0, are becoming 1 after taking mod 1
% px = mod(n - n_al * linspace(1,1,length(n)), 1);
px1 = mod(n1 - n_al * linspace(1,1,length(n1)), 1);

%For checking
% n1 - n_al
% mod(n1 - n_al, 1)
% graphics_toolkit('gnuplot');

%% Figure for old method
% figure(2);
% scatter(k_y/pi, px);
% hold on
% plot(k_y/pi, px, '--');
% xlabel('$k_y (\pi)$');
% ylabel('$p_x$ (modulo 1)');
% xlim([-1 1]);
% ylim([0 1]);
% title(strcat('$m = ',num2str(m_in_paper),'$'));
% grid on;

figure(1);
scatter(k_y/pi, px1);
hold on
plot(k_y/pi, px1, '--');
xlabel('$k_y (\pi)$', "fontsize", 30);
ylabel('$p_x$ (modulo 1)', "fontsize", 30);
xlim([-1 1]);
ylim([0 1]);
title(strcat('$m = ',num2str(m_in_paper),'$'), "fontsize", 30);
set(gca, "linewidth", 2, "fontsize", 30);
box on;

cd saved_plots/ssh_model
%define val_m_in_paper = '1pt75' etc.
filename = strcat('pol2DLx',num2str(m),'m',val_m_in_paper);

print(figure(1),'-dpdflatexstandalone', filename);
system(strcat("pdflatex\t", filename));
cd ../../

end

%%%1D%%%
function polarization1D(nsites)
clear all;
close all;
t = 1;
t1 = 1;
nsites = 100;
no_of_points = 30;
m_in_paper = linspace(-2,2,no_of_points);
m = nsites; n = 1; p = 1;
lattices;
angmom;
for i=1:no_of_points
	h_ssh = kron(m_in_paper(i) * M1D + t1 * CX1Dnp, sigma_x) + t * kron(SX1Dnp, sigma_y);
	[eigenstates, energy_eigenvalues] = eig(h_ssh);
	%value of x for 2k and 2k - 1 th eigenstate is k
	xmatrix = kron(diag(linspace(1,m,m)), eye(2));
	%reduced eigenstates
	filled_eigenstates = eigenstates(:, 1:m);
	x_exp = transpose(filled_eigenstates) * expm((2 * pi * j/m) * xmatrix) * conj(filled_eigenstates);
	n1(i) = real((-j/(2*pi)) * trace(logm(diag(eig(x_exp)))));
	%average value of number_at is (0 + 1)/2 = 1/2; Then n_al = trace(1/2 * xmatrix)/m = 1/2 * 2 * m * (m+1)/m = (m+1)/2;

end

n_al = (m+1)/2;
n1 = round(n1 * 100)/100; %This will round it of to 2 decimal places, otherwise values like -2e-5, 
n2 = round(n2 * 100)/100; %which should really be 0, are becoming 1 after taking mod 1
px1 = mod(n1 - n_al * linspace(1,1,length(n1)), 1);
px2 = mod(n2 - n_al * linspace(1,1,length(n2)), 1);

% graphics_toolkit('gnuplot');

figure 1;
scatter(m_in_paper,px1);
hold on
plot(m_in_paper, px1, '--');
xlabel('$m$');
ylabel('$p_x$ (modulo 1)');
grid on;

end

function Zakphase()
	m_in_paper = linspace(-2,2,100);
	k = linspace(0,2 * pi,20000);
	binsize = 2*pi/length(k);
 % Phase = (-i/2 pi) \int_{0}{2\pi} \frac{h'(k)}{h(k)} dk, where h(k) = m + t exp(ik), t is set to be 1
	for p = 1:length(m_in_paper)
		result(p) = 0;
		for q = 1:length(k)
			result(p) += real((1/(2*pi)) * (exp(j * k(q))/(m_in_paper(p) + exp(j * k(q)))) * binsize);
		end
	end
	figure(1)
	plot(m_in_paper, result);
	xlabel('$m$','fontsize',30);
	ylabel('Zak Phase','fontsize',30);
	set(gca, "linewidth", 2, "fontsize", 20);
	cd saved_plots/ssh_model/Zak
	print(figure(1), '-dpdflatexstandalone','Zakphase');
	system('pdflatex Zakphase');
	system('rm *.log *.aux');
	cd ../../..
end