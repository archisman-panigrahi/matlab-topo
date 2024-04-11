 t1 = 1;
function Energyspectra(Lx,Ly,npts)
close all;
tic;
m = Lx; n = Ly; p = 1;
lattices2D;
angmom;
m_in_paper = linspace (-4,4,npts);
t0 = 1;
t1 = 1;
for i = 1: npts
	h_cinp = t1 * (kron(SX2Dnp, sigma_x) + kron(SY2Dnp, sigma_y)) - t0 * kron(m_in_paper(i) * M2D - CX2Dnp - CY2Dnp, sigma_z);
	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) - t0 * kron(m_in_paper(i) * M2D - CX2Dp - CY2Dp, sigma_z);

	energy_eigenvalues_np = eig(h_cinp);
	energy_eigenvalues_p = eig(h_cip);
	for r=1:2*m*n
		energies_np( (i-1)*(2*m*n) + r ) = energy_eigenvalues_np(r);
		energies_p( (i-1)*(2*m*n) + r ) = energy_eigenvalues_p(r);
	end
end

xvar = kron(m_in_paper,linspace(1,1,2*m*n));
toc;
%graphics_toolkit('qt');
text("interpreter","latex");

figure(1);
scatter(xvar,energies_np,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for OBC", 'fontsize', 30)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

figure(2);
scatter(xvar,energies_p,'.');
xlabel('$m$', 'fontsize', 20);
ylabel('$E$', 'fontsize', 20);
title("Allowed energies for PBC", 'fontsize', 30)
set(gca, "linewidth", 2, "fontsize", 20);
axis tight;
box on;

% cd saved_plots/chern_insu/allowed_energies
% print(figure(1), '-dpdflatexstandalone', 'energiesOBC');
% print(figure(2), '-dpdflatexstandalone', 'energiesPBC');
% system("pdflatex energiesOBC");
% system("pdflatex energiesPBC");
% system("rm *.log *.aux")
% cd ../../..

end

function showEigenvalues(nsitex, nsitey, m_in_paper)
close all;
m = nsitex;
n = nsitey;
p = 1;
lattices2D;
angmom;
t0 = 1;
t1 = 1;

	h_cinp = t1 * (kron(SX2Dnp, sigma_x) + kron(SY2Dnp, sigma_y)) - t0 * kron(m_in_paper * M2D - CX2Dnp - CY2Dnp, sigma_z);
	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) - t0 * kron(m_in_paper * M2D - CX2Dp - CY2Dp, sigma_z);
	clear S*p C*p M2D;
	energy_eigenvalues_np = eig(h_cinp);
	energy_eigenvalues_p = eig(h_cip);
	clear h_cip h_cinp;
	%graphics_toolkit('gnuplot');
	text("interpreter","latex");

	figure(1);
	scatter(linspace(1,2*m*n,2*m*n), energy_eigenvalues_np);
	ylabel("$E$", 'fontsize', 20);
	title("OBC", 'fontsize', 20);
	axis([0 2*m*n]);
	set(gca, "linewidth", 2, "fontsize", 20);
	box on;

	figure(2);
	scatter(linspace(1,2*m*n,2*m*n), energy_eigenvalues_p);
	ylabel("$E$", 'fontsize', 20);
	title("PBC", 'fontsize', 20);
	axis([0 2*m*n]);
	set(gca, "linewidth", 2, "fontsize", 20);
	box on;

	cd saved_plots/chern_insu/eigenvalues
	mkdir(num2str(m_in_paper))
	cd(num2str(m_in_paper))
	print(figure(1), '-dpdflatexstandalone', 'energiesOBC');
	print(figure(2), '-dpdflatexstandalone', 'energiesPBC');
	system("pdflatex energiesOBC");
	system("pdflatex energiesPBC");
	system("rm *.log *.aux")
	cd ../../../..

end

function showEdgestate(nsitex, nsitey, m_in_paper)
	close all;
	t0 = 1;
	t1 = 1;
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;
	h_cinp = t1 * (kron(SX2Dnp, sigma_x) + kron(SY2Dnp, sigma_y)) - t0 * kron(m_in_paper * M2D - CX2Dnp - CY2Dnp, sigma_z);
	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) - t0 * kron(m_in_paper * M2D - CX2Dp - CY2Dp, sigma_z);
	clear S*p C*p M2D; 
	[eigenstates_np,energy_eigenvalues_np] = eig(h_cinp);
	[eigenstates_p,energy_eigenvalues_p] = eig(h_cip);
	clear h_cip h_cinp;
	high_occ_np = eigenstates_np(:, m*n);
	high_occ_p = eigenstates_p(:, m*n);
	clear eigenstates_p eigenstates_np;

	for k = 1:m*n
		prob_density_np(k) = abs(high_occ_np(2*k - 1))^2 + abs(high_occ_np(2*k))^2;
		prob_density_p(k) = abs(high_occ_p(2*k - 1))^2 + abs(high_occ_p(2*k))^2; 
	end
	clear high_occ_p high_occ_np;

	for k = 1:m
		for l = 1:n
			prob_density_np_z(k,l) = prob_density_np((l-1)*m + k);
			prob_density_p_z(k,l) = prob_density_p((l-1)*m + k);
		end
	end
	[xx,yy] = meshgrid(linspace(1,m,m),linspace(1,n,n));
	figure(1);
	contourf(xx,yy,prob_density_np_z);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title('OBC','fontsize', 20);
	% zlabel('probability density','fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	figure(2);
	contourf(xx,yy,prob_density_p_z);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title('PBC','fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	cd saved_plots/chern_insu/highest_occ_state
	mkdir(strcat('lx',num2str(m),'ly',num2str(n),'m',num2str(m_in_paper)))
	cd(strcat('lx',num2str(m),'ly',num2str(n),'m',num2str(m_in_paper)))
	print(figure(1),'-dpdflatexstandalone','OBCcontour');
	print(figure(2),'-dpdflatexstandalone','PBCcontour');
	system("pdflatex OBCcontour");
	system("pdflatex PBCcontour");
	system("rm *.log *.aux");
	cd ../../../..
end

%%%%Bandstructure
%2D
function bandstructureCI()
	close all;
	t0 = 1;
	t1 = 1;
	% For this Hamiltonian X and Y are symmetric, so we do
	% Gamma(0,0) --> X(pi,0) --> M(pi,pi) --> Gamma(0,0)
	% h_cinp = t1 * (kron(SX2Dnp, sigma_x) + kron(SY2Dnp, sigma_y)) - t0 * kron(m_in_paper * M2D - CX2Dnp - CY2Dnp, sigma_z);

	x1 = linspace(0,pi,50);
	y1 = linspace(0,0,50);

	x2 = linspace(pi,pi,50);
	y2 = linspace(0,pi,50);

	x3 = linspace(pi,0,50);
	y3 = linspace(pi,0,50);

	k_x = [x1,x2,x3];
	k_y = [y1,y2,y3];

	m_in_paper_1 = 0.5 * linspace(1,1,length(k_x));
	m_in_paper_2 = -0.5 * linspace(1,1,length(k_x));

	m_in_paper_3 = 1 * linspace(1,1,length(k_x));
	m_in_paper_4 = -1 * linspace(1,1,length(k_x));

	m_in_paper_5 = 1.5 * linspace(1,1,length(k_x));
	m_in_paper_6 = -1.5 * linspace(1,1,length(k_x));

	m_in_paper_7 = 2 * linspace(1,1,length(k_x));
	m_in_paper_8 = -2 * linspace(1,1,length(k_x));

	m_in_paper_9 = 3 * linspace(1,1,length(k_x));
	m_in_paper_10 = -3 * linspace(1,1,length(k_x));

	m_in_paper_11 = 0 * linspace(1,1,length(k_x));

	E1 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_1 - cos(k_x) - cos(k_y))).^2 );
	E2 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_1 - cos(k_x) - cos(k_y))).^2 );
	E3 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_2 - cos(k_x) - cos(k_y))).^2 );
	E4 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_2 - cos(k_x) - cos(k_y))).^2 );

	E5 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_3 - cos(k_x) - cos(k_y))).^2 );
	E6 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_3 - cos(k_x) - cos(k_y))).^2 );
	E7 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_4 - cos(k_x) - cos(k_y))).^2 );
	E8 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_4 - cos(k_x) - cos(k_y))).^2 );

	E9 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_5 - cos(k_x) - cos(k_y))).^2 );
	E10 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_5 - cos(k_x) - cos(k_y))).^2 );
	E11 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_6 - cos(k_x) - cos(k_y))).^2 );
	E12 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_6 - cos(k_x) - cos(k_y))).^2 );

	E13 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_7 - cos(k_x) - cos(k_y))).^2 );
	E14 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_7 - cos(k_x) - cos(k_y))).^2 );
	E15 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_8 - cos(k_x) - cos(k_y))).^2 );
	E16 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_8 - cos(k_x) - cos(k_y))).^2 );

	E17 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_9 - cos(k_x) - cos(k_y))).^2 );
	E18 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_9 - cos(k_x) - cos(k_y))).^2 );
	E19 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_10 - cos(k_x) - cos(k_y))).^2 );
	E20 = - sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_10 - cos(k_x) - cos(k_y))).^2 );

	E21 = sqrt( (t1 * sin(k_x)).^2 + (t1 * sin(k_y)).^2 + (t0 * (m_in_paper_11 - cos(k_x) - cos(k_y))).^2 );

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
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-3.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-3.35, '$X$', 'fontsize', 30)
	text(100,-3.35, '$M$', 'fontsize', 30)
	text(150,-3.35, '$\Gamma$', 'fontsize', 30)

	figure(2);
	plot(linspace(1,length(k_x),length(k_x)), E5,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E6, 'linewidth', 4, 'color', 'blue');
	hold on
	plot(linspace(1,length(k_x),length(k_x)), E7,'linewidth', 4, 'color', 'red', linspace(1,length(k_x),length(k_x)), E8, 'linewidth', 4, 'color', 'red');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-3.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-3.35, '$X$', 'fontsize', 30)
	text(100,-3.35, '$M$', 'fontsize', 30)
	text(150,-3.35, '$\Gamma$', 'fontsize', 30)

	figure(3);
	plot(linspace(1,length(k_x),length(k_x)), E9,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E10, 'linewidth', 4, 'color', 'blue');
	hold on
	plot(linspace(1,length(k_x),length(k_x)), E11,'linewidth', 4, 'color', 'red', linspace(1,length(k_x),length(k_x)), E12, 'linewidth', 4, 'color', 'red');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-4.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-4.35, '$X$', 'fontsize', 30)
	text(100,-4.35, '$M$', 'fontsize', 30)
	text(150,-4.35, '$\Gamma$', 'fontsize', 30)

	figure(4);
	plot(linspace(1,length(k_x),length(k_x)), E13,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E14, 'linewidth', 4, 'color', 'blue');
	hold on
	plot(linspace(1,length(k_x),length(k_x)), E15,'linewidth', 4, 'color', 'red', linspace(1,length(k_x),length(k_x)), E16, 'linewidth', 4, 'color', 'red');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-4.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-4.35, '$X$', 'fontsize', 30)
	text(100,-4.35, '$M$', 'fontsize', 30)
	text(150,-4.35, '$\Gamma$', 'fontsize', 30)

	figure(5);
	plot(linspace(1,length(k_x),length(k_x)), E17,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E18, 'linewidth', 4, 'color', 'blue');
	hold on
	plot(linspace(1,length(k_x),length(k_x)), E19,'linewidth', 4, 'color', 'red', linspace(1,length(k_x),length(k_x)), E20, 'linewidth', 4, 'color', 'red');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-6.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-6.35, '$X$', 'fontsize', 30)
	text(100,-6.35, '$M$', 'fontsize', 30)
	text(150,-6.35, '$\Gamma$', 'fontsize', 30)

	figure(6);
	plot(linspace(1,length(k_x),length(k_x)), E21,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), -E21, 'linewidth', 4, 'color', 'blue');
	ylabel("$E$", 'fontsize', 30);
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 2, "fontsize", 20);

	text(1,-2.35, '$\Gamma$', 'fontsize', 30) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,-2.35, '$X$', 'fontsize', 30)
	text(100,-2.35, '$M$', 'fontsize', 30)
	text(150,-2.35, '$\Gamma$', 'fontsize', 30)
	%Save the plots
	cd saved_plots/chern_insu/bandstructure
	print(figure(1),"-dpdflatexstandalone","bandsm0pt5")
	print(figure(2),"-dpdflatexstandalone","bandsm1")
	print(figure(3),"-dpdflatexstandalone","bandsm1pt5")
	print(figure(4),"-dpdflatexstandalone","bandsm2")
	print(figure(5),"-dpdflatexstandalone","bandsm3")
	print(figure(6),"-dpdflatexstandalone","bandsm0")

	system("pdflatex bandsm0pt5")
	system("pdflatex bandsm1")
	system("pdflatex bandsm1pt5")
	system("pdflatex bandsm2")
	system("pdflatex bandsm3")
	system("pdflatex bandsm0")
	system("rm *.aux *.log")
	cd ../../..
end

function BottIndexOld(nsitex, nsitey)
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m_in_paper = linspace(-4,4,no_of_points);
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
		h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) - t0 * kron(m_in_paper(i) * M2D - CX2Dp - CY2Dp, sigma_z);

		[eigenstates, energy_eigenvalues] = eig(h_cip);
		clear h*p;

		filled_eigenstates = eigenstates(:,1:(m*n));
		clear eigenstates;

		bott = transpose(filled_eigenstates) * U_x * conj(filled_eigenstates) * transpose(filled_eigenstates) * V_y * conj(filled_eigenstates) * transpose(filled_eigenstates) * ctranspose(U_x) * conj(filled_eigenstates) * transpose(filled_eigenstates) * ctranspose(V_y) * conj(filled_eigenstates);
		n1(i) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end
 %This will round it of to 2 decimal places, otherwise values like -2e-5, 
 %which should really be 0, are becoming 1 after taking mod 1
 
% px1 = mod(round(n1 * 10000)/10000, 1);


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('$m$', "fontsize", 30);
	ylabel('Bott Index', "fontsize", 30);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	% figure(2);
	% scatter(m_in_paper, px1);
	% hold on
	% plot(m_in_paper, px1, '--');
	% xlabel('$m$', "fontsize", 30);
	% ylabel('Bott Index', "fontsize", 30);
	% set(gca, "linewidth", 2, "fontsize", 30);
	% box on;

	cd saved_plots/chern_insu/Bott_Index
	filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% filename2 = strcat('modbottlx',num2str(m),'ly',num2str(n));
	print(figure(1),'-dpdflatexstandalone', filename1);
	% print(figure(2),'-dpdflatexstandalone', filename2);
	system(strcat("pdflatex\t", filename1));
	% system(strcat("pdflatex\t", filename2));
	system("rm *.log *.aux")
	cd ../../..

end

function BottIndexShow(nsitex, nsitey, m_in_paper)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);

	[eigenstates, energy_eigenvalues] = eig(h_cip);
	system("rm hamiltonianReal.dat");
	system("rm hamiltonianImag.dat");
	dlmwrite('hamiltonianReal.dat', real(h_cip), ' ');
	dlmwrite('hamiltonianImag.dat', imag(h_cip), ' ');
	clear h_cip;

	filled_eigenstates = eigenstates(:,1:(m*n));
	clear eigenstates;
	%Without conj(), it gives negative of Chern number. Which one is the convention?
	P = filled_eigenstates * ctranspose(filled_eigenstates);
	dlmwrite('PReal.dat', real(P), ' ');
	dlmwrite('PImag.dat', imag(P), ' ');
	clear filled_eigenstates
	U = P * U_x * P + (eye(2*m*n) - P);
	V = P * V_y * P + (eye(2*m*n) - P);
	bott = U * V * ctranspose(U) * ctranspose(V);
	bottindex = real((j/(2*pi)) * trace(logm(diag(eig(bott)))))
toc;
end

function BottIndex(nsitex, nsitey)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);

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
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	cd saved_plots/chern_insu/Bott_Index
	filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	print(figure(1),'-dpdflatexstandalone', filename1);
	system(strcat("pdflatex\t", filename1));
	system("rm *.log *.aux")
	cd ../../..
	toc;
end

function BottIndexOrder2(nsitex, nsitey)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(CX2Dp - CY2Dp, sigma_x) + kron(SXpYp2D, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);

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
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	# cd saved_plots/chern_insu/Bott_Index
	# filename1 = strcat('bott2ndorderlx',num2str(m),'ly',num2str(n));
	# print(figure(1),'-dpdflatexstandalone', filename1);
	# system(strcat("pdflatex\t", filename1));
	# system("rm *.log *.aux")
	# cd ../../..
	# toc;
end

function BottIndexOrder3(nsitex, nsitey)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(3*SXpCYp2D - (S2X2Dp/2) - 2*SX2Dp, sigma_x) + kron(3*CXpSYp2D - (S2Y2Dp/2) - 2*SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);
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
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end

	toc;
	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	# cd saved_plots/chern_insu/Bott_Index
	# filename1 = strcat('bott3rdorderlx',num2str(m),'ly',num2str(n));
	# print(figure(1),'-dpdflatexstandalone', filename1);
	# system(strcat("pdflatex\t", filename1));
	# system("rm *.log *.aux")
	# cd ../../..
end

function Order3Energies(nsitex, nsitey, m_in_paper)
	t1 = 1;
	t0 = 1;
	m = nsitex;
	n = nsitey;
	lattices2D;
	angmom;
	h_cip = t1 * (kron(3*SXpCYp2D - (S2X2Dp/2) - 2*SX2Dp, sigma_x) + kron(3*CXpSYp2D - (S2Y2Dp/2) - 2*SY2Dp, sigma_y)) - t0 * kron(m_in_paper(a) * M2D - CX2Dp - CY2Dp, sigma_z);

	energies = eig(h_cip);
	xvar = 1:(2*m*n);
	scatter(xvar, energies)
end

function Hermitian_check(a)
	b = size(a)
	c = 0;
	for xi = 1:b(1)
		for yi = 1:b(2)
			if (a(xi,yi) - conj(a(yi,xi))) == 0
				c += 0;
			else
				c += 1;
				xi,yi
			end
		end
	end
	c
end

function showBottIndexOrder3(nsitex, nsitey, m_in_paper)
		tic;
	close all;
	t0 = 1;
	t1 = 1;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

		h_cip = t1 * (kron(3*SXpCYp2D - (S2X2Dp/2) - 2*SX2Dp, sigma_x) + kron(3*CXpSYp2D - (S2Y2Dp/2) - 2*SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);
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
		bott_index = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))))

	toc;
end

function BottIndexHermi(nsitex, nsitey, hx, hy, hz)
	tic;
	close all;
	t0 = 1;
	t1 = 1;
	hx,hy,hz
	no_of_points = 40;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron([1, 1],kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1])));
	ymatrix2D = diag(kron([1, 1],kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1])));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);
		h_Hermi = kron(sigma_x + j*sigma_y, h_cip)/2 + kron(sigma_x - j * sigma_y, ctranspose(h_cip))/2;
		[eigenstates, energy_eigenvalues] = eig(h_Hermi);
		clear h_cip h_Hermi;

		filled_eigenstates = eigenstates(:,1:(2*m*n));
		clear eigenstates;
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		clear filled_eigenstates
		U = P * U_x * P + (eye(4*m*n) - P);
		V = P * V_y * P + (eye(4*m*n) - P);
		bott = U * V * ctranspose(U) * ctranspose(V);
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	% cd saved_plots/chern_insu/Bott_Index
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	% toc;
end

function BottIndexHermiOrder2(nsitex, nsitey, hx, hy, hz)
	tic;
	close all;
	t0 = 1;
	t1 = 1;
	hx,hy,hz
	no_of_points = 40;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron([1, 1],kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1])));
	ymatrix2D = diag(kron([1, 1],kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1])));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(CX2Dp - CY2Dp, sigma_x) + kron(SXpYp2D, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);
		h_Hermi = kron(sigma_x + j*sigma_y, h_cip)/2 + kron(sigma_x - j * sigma_y, ctranspose(h_cip))/2;
		[eigenstates, energy_eigenvalues] = eig(h_Hermi);
		clear h_cip h_Hermi;

		filled_eigenstates = eigenstates(:,1:(2*m*n));
		clear eigenstates;
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		clear filled_eigenstates
		U = P * U_x * P + (eye(4*m*n) - P);
		V = P * V_y * P + (eye(4*m*n) - P);
		bott = U * V * ctranspose(U) * ctranspose(V);
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	% cd saved_plots/chern_insu/Bott_Index
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	% toc;
end

function BottIndexHermiOrder3(nsitex, nsitey, hx, hy, hz)
	tic;
	close all;
	t0 = 1;
	t1 = 1;
	hx,hy,hz
	no_of_points = 40;

	m_in_paper = linspace(-4,4,no_of_points);
	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron([1, 1],kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1])));
	ymatrix2D = diag(kron([1, 1],kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1])));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	for a=1:no_of_points
		h_cip = t1 * (kron(3*SXpCYp2D - (S2X2Dp/2) - 2*SX2Dp, sigma_x) + kron(3*CXpSYp2D - (S2Y2Dp/2) - 2*SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper(a) * M2D, sigma_z);
		h_Hermi = kron(sigma_x + j*sigma_y, h_cip)/2 + kron(sigma_x - j * sigma_y, ctranspose(h_cip))/2;
		[eigenstates, energy_eigenvalues] = eig(h_Hermi);
		clear h_cip h_Hermi;

		filled_eigenstates = eigenstates(:,1:(2*m*n));
		clear eigenstates;
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		clear filled_eigenstates
		U = P * U_x * P + (eye(4*m*n) - P);
		V = P * V_y * P + (eye(4*m*n) - P);
		bott = U * V * ctranspose(U) * ctranspose(V);
		n1(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bott)))));
	end


	figure(1);
	scatter(m_in_paper, n1);
	hold on
	plot(m_in_paper, n1, '--');
	xlabel('\textbf{m}', "fontsize", 40);
	ylabel('\textbf{Bott Index}', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 40);
	box on;

	% cd saved_plots/chern_insu/Bott_Index
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	% toc;
end

function [IPR_reg,IPR_special] = inverse_partitcipation_ratio(psi, nsitex, nsitey)
	IPR = 0;
	density = sum(reshape(abs(psi.^2), 2, []));
	IPR_reg=sum(density.^2);
	IPR_special=sum(density.^4)/sum(density.^2); %%Suggested by my friend Aaradhya Pandey
end

function BottIndexVsChemicalPotentialRandomOnsitePeriodic(nsitex, nsitey, m_in_paper, disorderpotential)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);

	random_array = 2*rand(m*n,1)-ones(m*n,1);
	zero_trace_random_array = random_array - (sum(random_array)/(m*n))*ones(m*n,1);

	disorderHamiltonian = disorderpotential * kron(diag(zero_trace_random_array), eye(2));
	%disorderHamiltonian = disorderpotential * kron(diag(2*rand(m*n,1))-eye(m*n), eye(2));
	disp(strcat("trace of disorderHamiltonian is ", num2str(trace(disorderHamiltonian))));
	[eigenstates, energy_eigenvalues] = eig(h_cip + disorderHamiltonian);

	for a=1:(2*m*n)
		clear filled_eigenstates P U V bottMatrix
		filled_eigenstates = eigenstates(:,1:a);
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bottMatrix = U * V * ctranspose(U) * ctranspose(V);
		BottIndex(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bottMatrix)))));
		disp(strcat("index = ",num2str(a), " mu = ",num2str(energy_eigenvalues(a,a)), " C = ",num2str(BottIndex(a)),"\n"))

		[IPR_array_reg(a),IPR_array_special(a)] = inverse_partitcipation_ratio(eigenstates(:,a), m, n);
	end


	figure(1);
	scatter(diag(energy_eigenvalues), BottIndex);
	hold on
	plot(diag(energy_eigenvalues), BottIndex, '--');
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('Bott Index', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(2);
	scatter(diag(energy_eigenvalues), IPR_array_reg);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR1', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(3);
	scatter(diag(energy_eigenvalues), IPR_array_special);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR2', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	% cd saved_plots/chern_insu/Bott_Index_disorder
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	toc;
end

function BottIndexOrder2VsChemicalPotentialRandomOnsitePeriodic(nsitex, nsitey, m_in_paper, disorderpotential)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(CX2Dp - CY2Dp, sigma_x) + kron(SXpYp2D, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);

	random_array = 2*rand(m*n,1)-ones(m*n,1);
	zero_trace_random_array = random_array - (sum(random_array)/(m*n))*ones(m*n,1);

	disorderHamiltonian = disorderpotential * kron(diag(zero_trace_random_array), eye(2));
	fprintf("Trace of disorder Hamiltonian is")
	trace(disorderHamiltonian)
	[eigenstates, energy_eigenvalues] = eig(h_cip + disorderHamiltonian);

	for a=1:(2*m*n)
		clear filled_eigenstates P U V bottMatrix
		filled_eigenstates = eigenstates(:,1:a);
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bottMatrix = U * V * ctranspose(U) * ctranspose(V);
		BottIndex(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bottMatrix)))));
		disp(strcat("index = ",num2str(a), " mu = ",num2str(energy_eigenvalues(a,a)), " C = ",num2str(BottIndex(a)),"\n"))

		[IPR_array_reg(a),IPR_array_special(a)] = inverse_partitcipation_ratio(eigenstates(:,a), m, n);
	end


	figure(1);
	scatter(diag(energy_eigenvalues), BottIndex);
	hold on
	plot(diag(energy_eigenvalues), BottIndex, '--');
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('$\mu$', "fontsize", 40);
	ylabel('Bott Index', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(2);
	scatter(diag(energy_eigenvalues), IPR_array_reg);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR1', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(3);
	scatter(diag(energy_eigenvalues), IPR_array_special);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR2', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	% cd saved_plots/chern_insu/Bott_Index_disorder
	% filename1 = strcat('Chern2','bottlx',num2str(m),'ly',num2str(n),'And100TimesM',num2str(round(100*m_in_paper)),'And100TimesdisorderAmp',num2str(round(100*disorderpotential)));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	toc;
end

function BottIndexOrder3VsChemicalPotentialRandomOnsitePeriodic(nsitex, nsitey, m_in_paper, disorderpotential)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(3*SXpCYp2D - (S2X2Dp/2) - 2*SX2Dp, sigma_x) + kron(3*CXpSYp2D - (S2Y2Dp/2) - 2*SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);
	random_array = 2*rand(m*n,1)-ones(m*n,1);
	zero_trace_random_array = random_array - (sum(random_array)/(m*n))*ones(m*n,1);

	disorderHamiltonian = disorderpotential * kron(diag(zero_trace_random_array), eye(2));
	fprintf("Trace of disorder Hamiltonian is")
	trace(disorderHamiltonian)
	%disorderHamiltonian = disorderpotential * kron(diag(2*rand(m*n,1))-eye(m*n), eye(2));
	[eigenstates, energy_eigenvalues] = eig(h_cip + disorderHamiltonian);

	for a=1:(2*m*n)
		clear filled_eigenstates P U V bottMatrix
		filled_eigenstates = eigenstates(:,1:a);
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bottMatrix = U * V * ctranspose(U) * ctranspose(V);
		BottIndex(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bottMatrix)))));
		disp(strcat("index = ",num2str(a), " mu = ",num2str(energy_eigenvalues(a,a)), " C = ",num2str(BottIndex(a)),"\n"))

		[IPR_array_reg(a),IPR_array_special(a)] = inverse_partitcipation_ratio(eigenstates(:,a), m, n);
	end


	figure(1);
	scatter(diag(energy_eigenvalues), BottIndex);
	hold on
	plot(diag(energy_eigenvalues), BottIndex, '--');
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('$\mu$', "fontsize", 40);
	ylabel('Bott Index', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(2);
	scatter(diag(energy_eigenvalues), IPR_array_reg);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR1', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	figure(3);
	scatter(diag(energy_eigenvalues), IPR_array_special);
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('IPR2', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	cd saved_plots/chern_insu/Bott_Index_disorder
	filename1 = strcat('Chern3','bottlx',num2str(m),'ly',num2str(n),'And100TimesM',num2str(round(100*m_in_paper)),'And100TimesdisorderAmp',num2str(round(100*disorderpotential)));
	print(figure(1),'-dpdflatexstandalone', filename1);
	system(strcat("pdflatex\t", filename1));
	system("rm *.log *.aux")
	cd ../../..
	toc;
end

function BottIndexVsChemicalPotentialRandomOnsiteOpen(nsitex, nsitey, m_in_paper, disorderpotential)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(SX2Dnp, sigma_x) + kron(SY2Dnp, sigma_y)) + t0 * kron(CX2Dnp + CY2Dnp - m_in_paper * M2D, sigma_z);
	disorderHamiltonian = disorderpotential * kron(diag(2*rand(m*n,1))-eye(m*n), eye(2));
	[eigenstates, energy_eigenvalues] = eig(h_cip + disorderHamiltonian);

	for a=1:(2*m*n)
		clear filled_eigenstates P U V bottMatrix
		filled_eigenstates = eigenstates(:,1:a);
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bottMatrix = U * V * ctranspose(U) * ctranspose(V);
		BottIndex(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bottMatrix)))));
		disp(strcat("index = ",num2str(a), " mu = ",num2str(energy_eigenvalues(a,a)), " C = ",num2str(BottIndex(a)),"\n"))
	end


	figure(1);
	scatter(diag(energy_eigenvalues), BottIndex);
	hold on
	plot(diag(energy_eigenvalues), BottIndex, '--');
	title(strcat("Disorder Potential w/t = ", num2str(disorderpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('Bott Index', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	% cd saved_plots/chern_insu/Bott_Index_disorder
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	toc;
end

function BottIndexVsChemicalPotentialQuasiPeriodic(nsitex, nsitey, m_in_paper, QPpotential, quasi)
	tic;
	close all;
	t0 = 1;
	t1 = 1;

	no_of_points = 80;

	m = nsitex;
	n = nsitey;
	p = 1;
	lattices2D;
	angmom;

	xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]));
	ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]));

	U_x = expm(2 * pi * j * xmatrix2D/m);
	V_y = expm(2 * pi * j * ymatrix2D/n);

	h_cip = t1 * (kron(SX2Dp, sigma_x) + kron(SY2Dp, sigma_y)) + t0 * kron(CX2Dp + CY2Dp - m_in_paper * M2D, sigma_z);
	QPMatrix = QPpotential * kron(diag(kron(linspace(1,1,n), sin(quasi*(1:m)))), eye(2));
	[eigenstates, energy_eigenvalues] = eig(h_cip + QPMatrix);

	for a=1:(2*m*n)
		clear filled_eigenstates P U V bottMatrix
		filled_eigenstates = eigenstates(:,1:a);
		%Without conj(), it gives negative of Chern number. Which one is the convention?
		P = conj(filled_eigenstates) * transpose(filled_eigenstates);
		U = P * U_x * P + (eye(2*m*n) - P);
		V = P * V_y * P + (eye(2*m*n) - P);
		bottMatrix = U * V * ctranspose(U) * ctranspose(V);
		BottIndex(a) = real((-j/(2*pi)) * trace(logm(diag(eig(bottMatrix)))));
		disp(strcat("index = ",num2str(a), " mu = ",num2str(energy_eigenvalues(a,a)), " C = ",num2str(BottIndex(a)),"\n"))
	end


	figure(1);
	scatter(diag(energy_eigenvalues), BottIndex);
	hold on
	plot(diag(energy_eigenvalues), BottIndex, '--');
	title(strcat("Quasiperiodic Potential w/t = ", num2str(QPpotential/t1)), "fontsize", 10)
	xlabel('\mu', "fontsize", 40);
	ylabel('Bott Index', "fontsize", 40);
	set(gca, "linewidth", 2, "fontsize", 30);
	box on;

	% cd saved_plots/chern_insu/Bott_Index_disorder
	% filename1 = strcat('bottlx',num2str(m),'ly',num2str(n));
	% print(figure(1),'-dpdflatexstandalone', filename1);
	% system(strcat("pdflatex\t", filename1));
	% system("rm *.log *.aux")
	% cd ../../..
	toc;
end

%%%%% Rewrite band structure with Non-Hermitian Terms
function NHBand_Structure(m_in_paper, hx, hy, hz,name1)
	close all;
	% For this Hamiltonian X and Y are symmetric, so we do
	% Gamma(0,0) --> X(pi,0) --> M(pi,pi) --> Gamma(0,0)

	x1 = linspace(0,pi,50);
	y1 = linspace(0,0,50);

	x2 = linspace(pi,pi,50);
	y2 = linspace(0,pi,50);

	x3 = linspace(pi,0,50);
	y3 = linspace(pi,0,50);

	k_x = [x1,x2,x3];
	k_y = [y1,y2,y3];

	m_in_paper_1 = m_in_paper * linspace(1,1,length(k_x));
	E1 = sqrt( (cos(k_x) + cos(k_y) - m_in_paper_1 + j*hz).^2 + (sin(k_x) +  j*hx).^2 + (sin(k_y) + j*hy).^2);
	E2 = sqrt( (cos(k_x) + cos(k_y) + m_in_paper_1 + j*hz).^2 + (sin(k_x) +  j*hx).^2 + (sin(k_y) + j*hy).^2);
	if((hx^2 + hy^2 + hz^2) == 0)
		E1 = real(E1);
		E2 = real(E2);
	else
		E1 = imag(E1);
		E2 = imag(E2);
	end

	figure(1);
	if (m_in_paper != 0)
		plot(linspace(1,length(k_x),length(k_x)), E1,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), E2, 'linewidth', 4, 'color', 'red');
		if((hx^2 + hy^2 + hz^2) == 0)
			hold on
			plot(linspace(1,length(k_x),length(k_x)), -E1,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), -E2, 'linewidth', 4, 'color', 'red');
		end
	else
		if((hx^2 + hy^2 + hz^2) == 0)
			plot(linspace(1,length(k_x),length(k_x)), E1,'linewidth', 4, 'color', 'blue', linspace(1,length(k_x),length(k_x)), -E1, 'linewidth', 4, 'color', 'blue');
		else
			plot(linspace(1,length(k_x),length(k_x)), E1,'linewidth', 4, 'color', 'blue');
		end
	end
	
	if (m_in_paper != 0)
		h = legend(strcat('$m=',num2str(m_in_paper),'$'),strcat('$m=',num2str(-m_in_paper),'$'));
	else
		h = legend(strcat('$m=',num2str(m_in_paper),'$'));
	end

	if((hx^2 + hy^2 + hz^2) == 0)
		ylabel("$E$", 'fontsize', 50);
		legend (h, "location", "north");
	else
		ylabel("$Im(E)$", 'fontsize', 50);
		legend (h, "location", "northeast");
	end
	set (h, "fontsize", 30, "color","none");
	% xlabel("$k_y(\pi)$)
	set(gca, "xgrid", "on");
	axis([1 length(k_x)]);
	set(gca,'xtick',[1 50 100 150]);
	set(gca, 'xticklabel',({'','','',''}));
	set(gca, "linewidth", 4, "fontsize", 50);
	axis tight;
	a = axis;
	where_to_put_labels = a(3) - 0.1 * (a(4) - a(3));
	text(1,where_to_put_labels, '$\Gamma$', 'fontsize', 40) %Note "$\Gamma$" won't work, and xticklabel cannot recognize \Gamma
	text(50,where_to_put_labels, '$X$', 'fontsize', 40)
	text(100,where_to_put_labels, '$M$', 'fontsize', 40)
	text(150,where_to_put_labels, '$\Gamma$', 'fontsize', 40)
	if((hx^2 + hy^2 + hz^2) == 0)
		text(110,((a(3) + a(4))/2)+((a(4) - a(3))/2)*0.2,strcat('$h_x=',num2str(hx),'$'),'fontsize',20)
		text(110,((a(3) + a(4))/2),strcat('$h_y=',num2str(hy),'$'),'fontsize',20)
		text(110,((a(3) + a(4))/2)-((a(4) - a(3))/2)*0.2,strcat('$h_z=',num2str(hz),'$'),'fontsize',20)
	else
		text(110,a(3)+((a(4) - a(3))/2)*0.6,strcat('$h_x=',num2str(hx),'$'),'fontsize',20)
		text(110,a(3)+((a(4) - a(3))/2)*0.4,strcat('$h_y=',num2str(hy),'$'),'fontsize',20)
		text(110,a(3)+((a(4) - a(3))/2)*0.2,strcat('$h_z=',num2str(hz),'$'),'fontsize',20)
	end
	cd saved_plots/chern_insu/
	mkdir('bandstructureNH');
	cd('bandstructureNH');
	print(figure(1),"-dpdflatexstandalone",name1);

	system(strcat('pdflatex',"\t",name1));
	system("rm *.log *.aux");
	cd ../../..
end

%%%%%%
% Semi infinite lattice
% Infinite in x direction, so that kx is a good quantum number, and open system in y
function kxGoodEnergySpectra(nsitey, m_in_paper)
	tic;
	close all;
	t0 = 1;
	t1 = 1;
	%hx,hy,hz
	no_of_points = 200;

	m = nsitey;
	lattices1D;
	angmom;

	kx = linspace(-pi,pi, no_of_points);
	for k = 1:length(kx)
		h_p = t0 * (sin(kx(k)) * kron(M1D, sigma_x) + kron(SX1Dp, sigma_y)) + t1 * kron( (m_in_paper - cos(kx(k))) * M1D - CX1Dp, sigma_z);
		h_np = t0 * (sin(kx(k)) * kron(M1D, sigma_x) + kron(SX1Dnp, sigma_y)) + t1 * kron( (m_in_paper - cos(kx(k))) * M1D - CX1Dnp, sigma_z);
 		
 		energies_p = eig(h_p);
 		energies_np = eig(h_np);

 		for l = 1:length(h_p)
 			A_p(l + (k-1)*length(h_p) ,1) = kx(k);
 			A_p(l + (k-1)*length(h_p) ,2) = energies_p(l);

 			A_np(l + (k-1)*length(h_np) ,1) = kx(k);
 			A_np(l + (k-1)*length(h_np) ,2) = energies_np(l);
 		end
 	end
 	figure()
 	scatter(A_p(:,1), A_p(:,2),'.')
 	xlabel('$k_x$','fontsize',20)
 	title('PBC Energy spectrum','fontsize',30)
	set(gca, "linewidth", 2, "fontsize", 20);
 	axis tight;

 	figure()
 	scatter(A_np(:,1), A_np(:,2),'.')
 	xlabel('$k_x$','fontsize',20)
 	title('OBC Energy spectrum','fontsize',30)
	set(gca, "linewidth", 2, "fontsize", 20);
 	axis tight;

	toc;
end

function ChernNumber(m, num, p, q)
	%Create a grid for kx, ky
	%Calculate the Eigenstates everywhere
	%Calculate Berry Curvatures
	%Add them
	tic;
	angmom;
	m, num
	kx = linspace(-pi,pi,num);
	ky = linspace(-pi,pi,num);

	deltaK = (2 * pi)/(num - 1)

	chernnumber1 = 0;
	chernnumber2 = 0;

	right_eigenvectors_1 = zeros(2,num * num);
	right_eigenvectors_2 = zeros(2,num * num);
	for a = 1:num
		for b = 1:num
			h = sin(kx(a)) * sigma_x +  sin(ky(b)) * sigma_y + (m - cos(kx(a)) - cos(ky(b))) * sigma_z;

			[right_eigvec, eigval, left_eigvec] = eig(h);

			right_eigenvectors_1(1,a + num * (b - 1)) = right_eigvec(1,1);
			right_eigenvectors_1(2,a + num * (b - 1)) = right_eigvec(2,1);

			right_eigenvectors_2(1,a + num * (b - 1)) = left_eigvec(1,2);
			right_eigenvectors_2(2,a + num * (b - 1)) = left_eigvec(2,2);
		end
	end

	% for a = 2:(num - 1)
	% 	for b = 2:(num - 1)
	% 		D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(:,(a-1) + num * (b - 1)))/(2 * deltaK);
	% 		D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b)) - right_eigenvectors_1(1,a + num * (b - 2)))/(2 * deltaK);

	% 		chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;
		
	% 		D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,(a-1) + num * (b - 1)))/(2 * deltaK);
	% 		D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b)) - right_eigenvectors_2(1,a + num * (b - 2)))/(2 * deltaK);

	% 		chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
	% 	end
	% end

	for a = 2:(num - 1)
		for b = 2:(num - 1)
			D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(:,a+ num * (b - 1)))/(deltaK);
			D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b)) - right_eigenvectors_1(:,a + num * (b - 1)))/(deltaK);

			if( a==p)
				if (b == q)
					this = right_eigenvectors_1(:,a+ num * (b - 1))
					right_one = right_eigenvectors_1(:,a + 1+ num * (b-1))
					upper_one = right_eigenvectors_1(:,a+ num * b)
					D_psi_1_kx
					D_psi_1_ky
					BerryCurvature = j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx)
				end
			end


			chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;
		
			D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(:,(a) + num * (b - 1)))/(deltaK);
			D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b)) - right_eigenvectors_2(:,a + num * (b - 1)))/(deltaK);

			chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
		end
	end

	% for a = 1:num
	% 	for b = 1:num
	% 		a, b
	% 		if (a==num)
	% 			D_psi_1_kx = (right_eigenvectors_1(:,1 + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 		else
	% 			D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 		end

	% 		if (b==num)
	% 			D_psi_1_ky = (right_eigenvectors_1(:,a) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 		else
	% 			D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 		end

	% 		chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;
		
	% 		if (a==num)
	% 			D_psi_2_kx = (right_eigenvectors_2(:,1 + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 		else
	% 			D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 		end

	% 		if (b==num)
	% 			D_psi_2_ky = (right_eigenvectors_2(:,a) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 		else
	% 			D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 		end
	% 		chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
	% 	end
	% end

	% for a = 1:num
	% 	for b = 1:num
	% 		a, b

	% 		D_psi_1_kx = (right_eigenvectors_1(:,mod(a+1, num) + num * mod((b - 1), num)) - right_eigenvectors_1(1,mod(a,num) + num * mod((b - 1),num)))/deltaK;
	% 		if(a == num)
	% 			D_psi_1_kx = (right_eigenvectors_1(:,1 + num * mod((b - 1), num)) - right_eigenvectors_1(1,mod(a,num) + num * mod((b - 1),num)))/deltaK;
	% 		end

	% 		D_psi_1_ky = (right_eigenvectors_1(:,mod(a, num) + num * (mod(b, num))) - right_eigenvectors_1(1,mod(a,num) + num * mod((b - 1),num)))/deltaK;

	% 		if (b == num)
	% 			D_psi_1_ky = (right_eigenvectors_1(:,mod(a, num) + num * (mod(b, num))) - right_eigenvectors_1(1,mod(a,num) + num * mod((b - 1),num)))/deltaK;

	% 		chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;
		
	% 		D_psi_2_kx = (right_eigenvectors_2(:,mod((a+1), num) + num * mod((b - 1), num)) - right_eigenvectors_2(1,mod(a, num) + num * mod((b - 1), num)))/deltaK;
	% 		D_psi_2_ky = (right_eigenvectors_2(:,mod(a, num) + num * (mod(b, num))) - right_eigenvectors_2(1,mod(a, num) + num * mod((b - 1), num)))/deltaK;

	% 		chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
	% 	end
	% end

	% for a = 1: num
	% 	for b = 1:num
	% 		D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(:,(a-1) + num * (b - 1)))/(2 * deltaK);
	% 		D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b)) - right_eigenvectors_1(1,a + num * (b - 2)))/(2 * deltaK);

	% 		chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;
		
	% 		D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,(a-1) + num * (b - 1)))/(2 * deltaK);
	% 		D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b)) - right_eigenvectors_2(1,a + num * (b - 2)))/(2 * deltaK);

	% 		chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
	% 	end
	% end

	% % horizontal lines, except right corners
	% for a = 1:(num - 1)
	% 	% lower horizontal line
	% 	b = 1;
	% 	D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% 	D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;

	% 	% upper horizontal line
	% 	b = num;
	% 	D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_1_ky = (right_eigenvectors_1(:,a) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% 	D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_2_ky = (right_eigenvectors_2(:,a) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;

	% end

	
	% for b = 2:(num - 1)
	% 	%left vertical line
	% 	a = 1;
	% 	D_psi_1_kx = (right_eigenvectors_1(:,(a+1) + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_1_ky = (right_eigenvectors_1(:,a + num * (b-1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% 	D_psi_2_kx = (right_eigenvectors_2(:,(a+1) + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_2_ky = (right_eigenvectors_2(:,a + num * (b-1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;
	% 	%right vertical line, except right corners
	% 	a = num;
	% 	D_psi_1_kx = (right_eigenvectors_1(:,1 + num * (b - 1)) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_1_ky = (right_eigenvectors_1(:,a + num * b) - right_eigenvectors_1(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% 	D_psi_2_kx = (right_eigenvectors_2(:,1 + num * (b - 1)) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	D_psi_2_ky = (right_eigenvectors_2(:,a + num * b) - right_eigenvectors_2(1,a + num * (b - 1)))/deltaK;
	% 	chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;

	% end


	% % right lower corner

	% D_psi_1_kx = (right_eigenvectors_1(:,1) - right_eigenvectors_1(:,num))/deltaK;
	% D_psi_1_ky = (right_eigenvectors_1(:, 2 * num) - right_eigenvectors_1(:,num))/deltaK;
	% chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% D_psi_2_kx = (right_eigenvectors_2(:,1) - right_eigenvectors_2(:,num))/deltaK;
	% D_psi_2_ky = (right_eigenvectors_2(:, 2 * num) - right_eigenvectors_2(:,num))/deltaK;
	% chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;

	% % right upper corner
	% D_psi_1_kx = (right_eigenvectors_1(:,num * (num - 1) + 1) - right_eigenvectors_1(:,num * num))/deltaK;
	% D_psi_1_ky = (right_eigenvectors_1(:,num) - right_eigenvectors_1(:,num * num))/deltaK;
	% chernnumber1 += j * (D_psi_1_kx' * D_psi_1_ky - D_psi_1_ky' * D_psi_1_kx) * (1/(2*pi)) * deltaK * deltaK;

	% D_psi_2_kx = (right_eigenvectors_2(:,num * (num - 1) + 1) - right_eigenvectors_2(:,num * num))/deltaK;
	% D_psi_2_ky = (right_eigenvectors_2(:,num) - right_eigenvectors_2(:,num * num))/deltaK;
	% chernnumber2 += j * (D_psi_2_kx' * D_psi_2_ky - D_psi_2_ky' * D_psi_2_kx) * (1/(2*pi)) * deltaK * deltaK;

	chernnumber1
	chernnumber2

	toc;
end