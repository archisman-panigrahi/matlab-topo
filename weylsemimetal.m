 t1 = 1;
function Energyspectra(nsitex, nsitey, nsitez)
close all;
m = nsitex; n = nsitey; p = nsitez;
lattices3D;
angmom;
m_in_paper = linspace (-2,6,150);
t0 = 1;
t = 1;
tz = 1;
clear SZ3Dnp SZ3Dp;
for i = 1:150
	h_wsmnp = t * (kron(SX3Dnp, sigma_x) + kron(SY3Dnp, sigma_y)) + kron((tz * CZ3Dnp - m_in_paper(i) * M3D + t0 * (2 * M3D - CX3Dnp - CY3Dnp)), sigma_z);
	h_wsmp = t * (kron(SX3Dp, sigma_x) + kron(SY3Dp, sigma_y)) + kron((tz * CZ3Dp - m_in_paper(i) * M3D + t0 * (2 * M3D - CX3Dp - CY3Dp)), sigma_z);

	energy_eigenvalues_np = eig(h_wsmnp);
	energy_eigenvalues_p = eig(h_wsmp);
	for r=1:2*m*n*p
		energies_np( (i-1)*(2*m*n*p) + r ) = energy_eigenvalues_np(r);
		energies_p( (i-1)*(2*m*n*p) + r ) = energy_eigenvalues_p(r);		
	end
end

xvar = kron(m_in_paper,linspace(1,1,2*m*n*p));
graphics_toolkit('gnuplot');
text("interpreter","latex");

figure(1);
scatter(xvar,energies_np,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for OBC", 'fontsize', 20)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

figure(2);
scatter(xvar,energies_p,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for PBC", 'fontsize', 20)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

cd saved_plots/WSM/allowed_energies
print(figure(1), '-dpdflatexstandalone', 'energiesOBC');
print(figure(2), '-dpdflatexstandalone', 'energiesPBC');
system("pdflatex energiesOBC");
system("pdflatex energiesPBC");
system("rm *.log *.aux")
cd ../../..

end 

function listEnergies(lx,ly,lz,m_in_paper)
close all;
m = lx; n = ly; p = lz;
lattices3D;
angmom;
t0 = 1;
t = 1;
tz = 0.5;
clear SZ3Dnp SZ3Dp;
	h_wsmnp = t * (kron(SX3Dnp, sigma_x) + kron(SY3Dnp, sigma_y)) + kron((tz * CZ3Dnp - m_in_paper * M3D + t0 * (2 * M3D - CX3Dnp - CY3Dnp)), sigma_z);
	h_wsmp = t * (kron(SX3Dp, sigma_x) + kron(SY3Dp, sigma_y)) + kron((tz * CZ3Dp - m_in_paper * M3D + t0 * (2 * M3D - CX3Dp - CY3Dp)), sigma_z);

	energy_eigenvalues_np = eig(h_wsmnp);
	energy_eigenvalues_p = eig(h_wsmp);

xvar = linspace(1,2*m*n*p,2*m*n*p);
graphics_toolkit('gnuplot');
text("interpreter","latex");

figure(1);
scatter(xvar,energy_eigenvalues_np,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for OBC", 'fontsize', 20)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

figure(2);
scatter(xvar,energy_eigenvalues_p,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for PBC", 'fontsize', 20)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

cd saved_plots/WSM/listEnergies
print(figure(1), '-dpdflatexstandalone', 'energiesOBC');
print(figure(2), '-dpdflatexstandalone', 'energiesPBC');
system("pdflatex energiesOBC");
system("pdflatex energiesPBC");
system("rm *.log *.aux")
cd ../../..

end 

function showEigenstate(lx,ly,lz,m_in_paper)
	close all;
	m = lx; n = ly; p = lz;
	lattices3D;
	angmom;
	t0 = 1;
	t = 1;
	tz = 1;
	clear SZ3Dnp SZ3Dp;
	h_wsmnp = t * (kron(SX3Dnp, sigma_x) + kron(SY3Dnp, sigma_y)) + kron((tz * CZ3Dnp - m_in_paper * M3D + t0 * (2 * M3D - CX3Dnp - CY3Dnp)), sigma_z);
	h_wsmp = t * (kron(SX3Dp, sigma_x) + kron(SY3Dp, sigma_y)) + kron((tz * CZ3Dp - m_in_paper * M3D + t0 * (2 * M3D - CX3Dp - CY3Dp)), sigma_z);

	clear C*p S*p;
	[states_np,energy_eigenvalues_np] = eig(h_wsmnp);
	[states_p, energy_eigenvalues_p] = eig(h_wsmp);
	highest_state_np = states_np(:,m*n*p);
	highest_state_p = states_p(:,m*n*p);
	% In 3D, (i,k,l) is mapped to (l-1)*m*n + (k-1) * m + i, where l varies from 1 to p, k from 1 to n, i from 1 to m
	for a = 1:m
		for b = 1:n
			for c = 1:p
				% (a,b) intentionally swapped
				density_np(b,a,c) = abs(highest_state_np(2*((c-1)*m*n + (b-1)*m + a)))^2 + abs(highest_state_np(2*((c-1)*m*n + (b-1)*m + a) - 1))^2;
				density_p(b,a,c) = abs(highest_state_p(2*((c-1)*m*n + (b-1)*m + a)))^2 + abs(highest_state_p(2*((c-1)*m*n + (b-1)*m + a) - 1))^2;
			end
		end
	end

	[xx,yy,zz] = meshgrid(linspace(1,m,m),linspace(1,n,n),linspace(1,p,p));

	figure(3)
	scatter3(xx(:),yy(:),zz(:),1000*density_np(:),1000*density_np(:),"fill");
	colormap(rainbow);
	box on;
	colorbar('fontsize',20);
	set(gca, "linewidth", 2, "fontsize", 20);
	title('OBC', 'fontsize', 30);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	zlabel('$z$','fontsize', 20);

	figure(4)
	scatter3(xx(:),yy(:),zz(:),1000*density_p(:),1000*density_p(:),"fill");
	colormap(rainbow);
	box on;
	colorbar('fontsize',20);
	set(gca, "linewidth", 2, "fontsize", 20);
	title('PBC', 'fontsize', 30)
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	zlabel('$z$','fontsize', 20);

	cd saved_plots/WSM/density_plots
	mkdir(num2str(m_in_paper));
	cd(num2str(m_in_paper));
	name1 = strcat('OBC','m',num2str(m_in_paper),'lx',num2str(m),'ly',num2str(n),'lz',num2str(p));
	name2 = strcat('PBC','m',num2str(m_in_paper),'lx',num2str(m),'ly',num2str(n),'lz',num2str(p));

	print(figure(3),'-dpdflatexstandalone',name1);
	print(figure(4),'-dpdflatexstandalone',name2);

	system(strcat("pdflatex\t",name1));
	system(strcat("pdflatex\t",name2));
	system('rm *.log *.aux');
	
	cd ../../../..
end

function bandstructureWSM(mz,name)
	close all;
	t0 = 1;
	t = 1;
	tz = 0.5;
	%M(pi,pi,pi) --> Gamma(0,0,0) -> Z(0,0,pi) -> M -> R1(pi,pi,0) -> Gamma -> X(pi,0,0) -> R2(pi,0,pi) -> Gamma

	x1 = linspace(pi, 0, 50);
	y1 = linspace(pi, 0, 50);
	z1 = linspace(pi, 0, 50);

	x2 = linspace(0, 0, 50);
	y2 = linspace(0, 0, 50);
	z2 = linspace(0, pi, 50);

	x3 = linspace(0, pi, 50);
	y3 = linspace(0, pi, 50);
	z3 = linspace(pi, pi, 50);

	x4 = linspace(pi, pi, 50);
	y4 = linspace(pi, pi, 50);
	z4 = linspace(pi, 0, 50);

	x5 = linspace(pi, 0, 50);
	y5 = linspace(pi, 0, 50);
	z5 = linspace(0, 0, 50);

	x6 = linspace(0, pi, 50);
	y6 = linspace(0, 0, 50);
	z6 = linspace(0, 0, 50);

	x7 = linspace(pi, pi, 50);
	y7 = linspace(0, 0, 50);
	z7 = linspace(0, pi, 50);

	x8 = linspace(pi, 0, 50);
	y8 = linspace(0, 0, 50);
	z8 = linspace(pi, 0, 50);

	%h = t sin(kx) sigma_x + t sin(ky) sigma_y + sigma_z * (tz cos(kz) - mz + t0(2 - cos(kx) - cos(ky))
	k_x = [x1,x2,x3,x4,x5,x6,x7,x8];
	k_y = [y1,y2,y3,y4,y5,y6,y7,y8];
	k_z = [z1,z2,z3,z4,z5,z6,z7,z8];
	E1 = hypot((t * sin(k_x)), hypot ((t * sin(k_y)), (tz * cos(k_z) - mz + t0 * (2 - cos(k_x) - cos(k_y)))));
	E2 = -E1;

	% text("interpreter","latex");
	figure(1);
	plot(linspace(1,length(k_x),length(k_x)), E1,'linewidth', 5, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E2, 'linewidth', 5, 'color', 'blue');
	hold on
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150 200 250 300 350 400]);
	set(gca, 'xticklabel',({'','','','','','','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);
	a = ylim;
	b = a(1) - 0.4;
	text(1,b, '$M$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,b, '$\Gamma$', 'fontsize', 30)
	text(100,b, '$Z$', 'fontsize', 30)
	text(150,b, '$M$', 'fontsize', 30)
	text(200,b, '$R_1$', 'fontsize', 30)
	text(250,b, '$\Gamma$', 'fontsize', 30)
	text(300,b, '$X$', 'fontsize', 30)
	text(350,b, '$R_2$', 'fontsize', 30)
	text(400,b, '$\Gamma$', 'fontsize', 30)

	cd saved_plots/WSM/bandstructure
	print(figure(1),"-dpdflatexstandalone",strcat("bandsmz", name));

	system(strcat('pdflatex bandsmz',name));
	system("rm *.aux *.log")
	cd ../../..
end

function FermiArc(lx,mz,cutoff,name)
	close all;
	%h = t sin(kx) sigma_x + t sin(ky) sigma_y + sigma_z * (tz cos(kz) - mz + t0(2 - cos(kx) - cos(ky))
	tic;
	t = 1;
	t0 = 1;
	tz = 0.5;
	no_of_pts = 40;
	k_y = linspace(-pi,pi,no_of_pts);
	k_z = linspace(-pi,pi,no_of_pts);
	m = lx;
	n = 1;
	p = 1;
	lattices1D;
	angmom;
	for q = 1:no_of_pts
		for r = 1:no_of_pts
			h_wsm = kron(t * SX1Dnp, sigma_x) + kron(t * sin(k_y(q)) * M1D, sigma_y) + kron((tz * cos(k_z(r)) - mz + t0 * (2 - cos(k_y(q)))) * M1D - t0 * CX1Dnp, sigma_z);
			[states,energies] = eig(h_wsm);

			highest_state = states(:, m);

			for a= 1:m
				density(q,a,r) = (abs(highest_state(2*a)))^2 + (abs(highest_state(2*a - 1)))^2;
				%This is not density(a,q,r), although a denotes index in x
				%See https://stackoverflow.com/a/62394044/3128341
			end
		end
	end
	toc;
	density(abs(density)<cutoff) = 0;
	x = linspace(1,m,m);
	[xx,yy,zz] = meshgrid(x,k_y/pi,k_z/pi);
	figure(1)
	scatter3(xx(:),yy(:),zz(:),10 * density(:),density(:),"fill");
	colormap(rainbow);
	colorbar('fontsize',20);
	grid off;
	box off;
	xlabel('$x$','fontsize', 20);
	ylabel('$k_y (\pi)$','fontsize', 20);
	zlabel('$k_z (\pi)$','fontsize', 20);
	xlim([1 m]);
	ylim([-1 1]);
	zlim([-1 1]);
	hold on
	plot3([1,m],[-1,-1],[-1,-1],'linewidth',2,'k')
	plot3([m,m],[-1,1],[-1,-1],'linewidth',2,'k')
	plot3([m,1],[1,1],[-1,-1],'linewidth',2,'k')
	plot3([1,1],[1,-1],[-1,-1],'linewidth',2,'k')
	plot3([1,m],[-1,-1],[1,1],'linewidth',2,'k')
	plot3([m,m],[-1,1],[1,1],'linewidth',2,'k')
	plot3([m,1],[1,1],[1,1],'linewidth',2,'k')
	plot3([1,1],[1,-1],[1,1],'linewidth',2,'k')
	plot3([1,1],[-1,-1],[-1,1],'linewidth',2,'k')
	plot3([1,1],[1,1],[-1,1],'linewidth',2,'k')
	plot3([m,m],[-1,-1],[-1,1],'linewidth',2,'k')
	plot3([m,m],[1,1],[-1,1],'linewidth',2,'k')
	set(gca, "linewidth", 2, "fontsize", 20);
	% cd saved_plots/WSM/FermiArc
	% print(figure(1),'-dpdflatexstandalone',name);
	% system(strcat("pdflatex\t",name));
	% system('rm *.log *.aux');
	% cd ../../..
end

function FermiArc1Dvisualization(lx,mz,k_y,k_z)
	close all;
	t = 1;
	t0 = 1;
	tz = 0.5;
	m = lx;
	lattices1D;
	angmom;
	h_wsm = kron(t * SX1Dnp, sigma_x) + kron(t * sin(k_y) * M1D, sigma_y) + kron((tz * cos(k_z) - mz + t0 * (2 - cos(k_y))) * M1D - t0 * CX1Dnp, sigma_z);
	[states,energies] = eig(h_wsm);

	highest_state = states(:, m);
	length(highest_state)
	[maxval maxind] = max(abs(highest_state))
	for a= 1:m
		density(a) = (abs(highest_state(2*a)))^2 + (abs(highest_state(2*a - 1)))^2;
	end
	[maxval maxind] = max(abs(density))

	scatter(linspace(1,m,m),density);
end

function TopSurfacekykz(lx,mz,name)
	close all;
	%h = t sin(kx) sigma_x + t sin(ky) sigma_y + sigma_z * (tz cos(kz) - mz + t0(2 - cos(kx) - cos(ky))
	tic;
	t = 1;
	t0 = 1;
	tz = 0.5;
	no_of_pts = 100;
	%These needs to be accordingly changed to get detailed plots for different values of mz
	k_y = linspace(-pi,pi,no_of_pts);
	k_z = linspace(-pi,pi,no_of_pts);
	m = lx;
	lattices1D;
	angmom;
	for q = 1:length(k_y)
		for r = 1:length(k_z)
			h_wsm = kron(t * SX1Dnp, sigma_x) + kron(t * sin(k_y(q)) * M1D, sigma_y) + kron((tz * cos(k_z(r)) - mz + t0 * (2 - cos(k_y(q)))) * M1D - t0 * CX1Dnp, sigma_z);
			[states,energies] = eig(h_wsm);

			highest_state = states(:, m);
			density(q,r) = abs(highest_state(2*m))^2 + abs(highest_state(2*m - 1))^2 + abs(highest_state(1))^2 + abs(highest_state(2))^2;
			energy(q,r) = energies(m,m);
		end
	end
	max_density = max(max(density))
	[a,b] = find(density == max_density);
	nrg = energy(a,b)(1,1)
	toc;
	tic;
	[xx yy] = meshgrid(k_z,k_y);
	figure(1)
	contourf(xx/pi,yy/pi,density/max_density)
	colorbar('fontsize',20);
	xlabel('$k_z (\pi)$','fontsize',15) %Otherwise it does not fit inside the plot
	ylabel('$k_y (\pi)$','fontsize',20)
	title(strcat('$E = ',num2str(nrg),'$'),'fontsize',20)
	set(gca, "linewidth", 2, "fontsize", 20);
	cd saved_plots/WSM/topsurface
	print(figure(1),'-dpdflatexstandalone',name);
	system(strcat("pdflatex\t",name));
	system('rm *.log *.aux');
	cd ../../..
	toc;
end
function BottIndex(nsitex, nsitey, mz)
	close all;
	t0 = 1;
	t = 1;
	tz = 0.5;
	no_of_points = 80;
	k_z = linspace(-pi,pi,no_of_points);
	m_in_paper = mz;
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for i=1:no_of_points
		h_cip = t * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) - kron((tz * cos(k_z(i)) - m_in_paper) * M2D + t0 * 2 * M2D - t0 * (CX2Dp + CY2Dp), sigma_z);

		[eigenstates, energy_eigenvalues] = eig(h_cip);
		clear h*p;

		filled_eigenstates = eigenstates(:,1:(m*n));
		clear eigenstates;
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		clear filled_eigenstates
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bott = U * V * ctranspose(U) * ctranspose(V);
		n1(i) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end
	n1 = round(n1 * 100)/100;

	figure(2);
	scatter(k_z/pi, n1);
	hold on
	plot(k_z/pi, n1, '--');
	xlabel('$k_z (\pi)$', "fontsize", 30);
	ylabel('Bott Index', "fontsize", 30);
	title(strcat('$m_z = $',num2str(m_in_paper)),'fontsize',30);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;
	axis auto;
	cd saved_plots/WSM/Bott_Index
	filename1 = strcat('bottlx',num2str(m),'ly',num2str(n),'mz',num2str(m_in_paper));
	% filename1 = 'bottlx10ly10mz-1pt5' %'bottlx10ly10mz5pt5'
	print(figure(2),'-dpdflatexstandalone', filename1);
	system(strcat("pdflatex	\t", filename1));
	system("rm *.log *.aux");
	cd ../../..

end


function PrintBottIndex()

	mz = [0,1,2,3,4]; %Do mz = -1.5 and 5.5 manually because if filename has "." system will show error
	for a = 1:5
		BottIndex(10,10,mz(a));
	end
end

function Energy_in_BandStructure(k_x,k_y,k_z,m_z,t,t0,tz)
	k_x,k_y,k_z,t,t0,tz,m_z

	% E1 = hypot((t * sin(k_x)), hypot ((t * sin(k_y)), (tz * cos(k_z) - m_z + t0 * (2 - cos(k_y) - cos(k_z)))))
	E1 = t * sin(k_x)
	E2 = t * sin(k_y)
	E3 = tz * cos(k_z) - m_z + t0 * (2 - cos(k_x) - cos(k_y))
	E = hypot(E1, E2, E3)
end