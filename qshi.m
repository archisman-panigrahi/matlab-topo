	angmom;
	global gamma_1 = kron(sigma_z, sigma_x);
	global gamma_2 = kron(sigma_0, sigma_y);
	global gamma_3 = kron(sigma_0, sigma_z);
	global gamma_4 = kron(sigma_x, sigma_x);
	global gamma_5 = kron(sigma_y, sigma_x);

function EnergySpectrum_without_mass(nsitex, nsitey, m_in_paper)
	close all;
	global gamma_1;
	global gamma_2;
	global gamma_3;
	global gamma_4;
	global gamma_5;
	t1 = 1;
	t0 = 1;
	m = nsitex;
	n = nsitey;
	lattices2D;
	h_SHI_p = t1 * (kron(SX2Dp, gamma_1) + kron(SY2Dp, gamma_2)) - t0 * (kron(m_in_paper * M2D - CX2Dp - CY2Dp, gamma_3));
	h_SHI_np = t1 * (kron(SX2Dnp, gamma_1) + kron(SY2Dnp, gamma_2)) - t0 * (kron(m_in_paper * M2D - CX2Dnp - CY2Dnp, gamma_3));
	clear SX2Dp SY2Dp SX2Dnp SY2Dnp M2D;
	[states_p, energies_p] = eig(h_SHI_p);
	[states_np, energies_np] = eig(h_SHI_np);

	%Since there are total 4*m*n states
	highest_state_p1 = states_p(:, 2*m*n);
	highest_state_p2 = states_p(:, 2*m*n - 1);

	highest_state_np1 = states_np(:, 2*m*n);
	highest_state_np2 = states_np(:, 2*m*n - 1);

	energies_p = diag(energies_p);
	energies_np = diag(energies_np);
	clear states_p states_np;

	for a = 1:m
		for b = 1:n
			%(Indices are inverted because meshgrid inverts them)
			density_np1(b,a) = abs(highest_state_np1(4*((b-1)*m + a)))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 1))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 2))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 3))^2;
			density_np2(b,a) = abs(highest_state_np2(4*((b-1)*m + a)))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 1))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 2))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 3))^2;

			density_p1(b,a) = abs(highest_state_p1(4*((b-1)*m + a)))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 1))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 2))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 3))^2;
			density_p2(b,a) = abs(highest_state_p2(4*((b-1)*m + a)))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 1))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 2))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 3))^2;
		end
	end

	density_np = density_np1 + density_np2;
	density_p = density_p1 + density_p2;

	[xx,yy] = meshgrid(linspace(1,m,m),linspace(1,n,n));

	xvar = linspace(1,length(energies_p),length(energies_p));

	% graphics_toolkit('gnuplot');
	% text("interpreter","latex");

	figure(1);
	scatter(xvar,energies_p,'.');
	xlabel('$n$', 'fontsize', 20);
	ylabel('$E$', 'fontsize', 20);
	title("Allowed energies for PBC", 'fontsize', 20)
	set(gca, "linewidth", 2, "fontsize", 20);
	axis tight;
	box on;

	figure(2);
	scatter(xvar,energies_np,'.');
	xlabel('$n$', 'fontsize', 20);
	ylabel('$E$', 'fontsize', 20);
	title("Allowed energies for OBC", 'fontsize', 20)
	set(gca, "linewidth", 2, "fontsize", 20);
	axis tight;
	box on;

	title3 = strcat('PBC,','$E = ',num2str(energies_p(2*m*n)),'$')
	title4 = strcat('OBC,','$E = ',num2str(energies_np(2*m*n)),'$')

	figure(3);
	contourf(xx,yy,density_p);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title(title3,'fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	figure(4);
	contourf(xx,yy,density_np);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title(title4,'fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	cd saved_plots/QSHI/allowed_energies_without_mass
	print(figure(1), '-dpdflatexstandalone', 'energiesPBC');
	print(figure(2), '-dpdflatexstandalone', 'energiesOBC');
	system("pdflatex energiesOBC");
	system("pdflatex energiesPBC");
	system("rm *.log *.aux")
	cd ../../..

	cd saved_plots/QSHI/density_without_mass
	print(figure(3), '-dpdflatexstandalone', 'densityPBC');
	print(figure(4), '-dpdflatexstandalone', 'densityOBC');
	system("pdflatex densityOBC");
	system("pdflatex densityPBC");
	system("rm *.log *.aux")
	cd ../../..
end

function EnergySpectrum_with_mass(nsitex, nsitey, m_in_paper)
	close all;
	global gamma_1;
	global gamma_2;
	global gamma_3;
	global gamma_4;
	global gamma_5;
	t1 = 1;
	t0 = 1;
	delta = 0.3;
	m = nsitex;
	n = nsitey;
	lattices2D;
	h_SHI_p = t1 * (kron(SX2Dp, gamma_1) + kron(SY2Dp, gamma_2)) - t0 * (kron(m_in_paper * M2D - CX2Dp - CY2Dp, gamma_3)) + delta * kron(CX2Dp - CY2Dp, gamma_4);
	h_SHI_np = t1 * (kron(SX2Dnp, gamma_1) + kron(SY2Dnp, gamma_2)) - t0 * (kron(m_in_paper * M2D - CX2Dnp - CY2Dnp, gamma_3)) + delta * kron(CX2Dnp - CY2Dnp, gamma_4);
	clear SX2Dp SY2Dp SX2Dnp SY2Dnp M2D;
	[states_p, energies_p] = eig(h_SHI_p);
	[states_np, energies_np] = eig(h_SHI_np);

	%Since there are total 4*m*n states
	highest_state_p1 = states_p(:, 2*m*n);
	highest_state_p2 = states_p(:, 2*m*n - 1);

	highest_state_np1 = states_np(:, 2*m*n);
	highest_state_np2 = states_np(:, 2*m*n - 1);

	energies_p = diag(energies_p);
	energies_np = diag(energies_np);
	clear states_p states_np;

	for a = 1:m
		for b = 1:n
			%(Indices are inverted because meshgrid inverts them)
			density_np1(b,a) = abs(highest_state_np1(4*((b-1)*m + a)))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 1))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 2))^2 + abs(highest_state_np1(4*((b-1)*m + a) - 3))^2;
			density_np2(b,a) = abs(highest_state_np2(4*((b-1)*m + a)))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 1))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 2))^2 + abs(highest_state_np2(4*((b-1)*m + a) - 3))^2;

			density_p1(b,a) = abs(highest_state_p1(4*((b-1)*m + a)))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 1))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 2))^2 + abs(highest_state_p1(4*((b-1)*m + a) - 3))^2;
			density_p2(b,a) = abs(highest_state_p2(4*((b-1)*m + a)))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 1))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 2))^2 + abs(highest_state_p2(4*((b-1)*m + a) - 3))^2;
		end
	end

	density_np = density_np1 + density_np2;
	density_p = density_p1 + density_p2;

	[xx,yy] = meshgrid(linspace(1,m,m),linspace(1,n,n));

	xvar = linspace(1,length(energies_p),length(energies_p));

	% graphics_toolkit('gnuplot');
	% text("interpreter","latex");

	figure(1);
	scatter(xvar,energies_p,'.');
	xlabel('$n$', 'fontsize', 20);
	ylabel('$E$', 'fontsize', 20);
	title("Allowed energies for PBC", 'fontsize', 20)
	set(gca, "linewidth", 2, "fontsize", 20);
	axis tight;
	box on;

	figure(2);
	scatter(xvar,energies_np,'.');
	xlabel('$n$', 'fontsize', 20);
	ylabel('$E$', 'fontsize', 20);
	title("Allowed energies for OBC", 'fontsize', 20)
	set(gca, "linewidth", 2, "fontsize", 20);
	axis tight;
	box on;

	title3 = strcat('PBC,','$E = ',num2str(energies_p(2*m*n)),'$')
	title4 = strcat('OBC,','$E = ',num2str(energies_np(2*m*n)),'$')

	figure(3);
	contourf(xx,yy,density_p);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title(title3,'fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	figure(4);
	contourf(xx,yy,density_np);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	title(title4,'fontsize', 20);
	colorbar('fontsize',20);
	axis equal;
	axis tight;
	box on;
	set(gca, "linewidth", 2, "fontsize", 20);

	cd saved_plots/QSHI/allowed_energies_with_mass
	print(figure(1), '-dpdflatexstandalone', 'energiesPBC');
	print(figure(2), '-dpdflatexstandalone', 'energiesOBC');
	system("pdflatex energiesOBC");
	system("pdflatex energiesPBC");
	system("rm *.log *.aux")
	cd ../../..

	cd saved_plots/QSHI/density_with_mass
	print(figure(3), '-dpdflatexstandalone', 'densityPBC');
	print(figure(4), '-dpdflatexstandalone', 'densityOBC');
	system("pdflatex densityOBC");
	system("pdflatex densityPBC");
	system("rm *.log *.aux")
	cd ../../..
end