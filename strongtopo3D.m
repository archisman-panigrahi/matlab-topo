angmom;
global gamma1 = kron(sigma_z,sigma_0);
global gamma2 = kron(sigma_x,sigma_x);
global gamma3 = kron(sigma_x,sigma_y);
global gamma4 = kron(sigma_x,sigma_z);
global gamma1_prime = kron(sigma_y,sigma_0);
function showStatesAndEnergies(lx,ly,lz,cutoff,topological_order)
	%e.g. showStatesAndEnergies(4,4,4,0.5,1)
	global gamma1;
	global gamma2;
	global gamma3;
	global gamma4;
	global gamma1_prime;
	close all;
	tic;
	q = 0.5;
	m = lx; n = ly; p = lz;
	lattices3Dnp;
	t = 1;
	eta = 1;
	M = -2;
	delta1 = 0.3;
	delta2 = 0.3;

	h_dsm_np = kron(M * M3D + t * (CX3Dnp + CY3Dnp + CZ3Dnp), gamma1) + kron(eta * SX3Dnp, gamma2) + kron(eta * SY3Dnp, gamma3) + kron(eta * SZ3Dnp, gamma4);
	if topological_order == 2
		h_dsm_np += delta1 * kron(CX3Dnp - CY3Dnp, gamma1_prime);
	% elseif topological_order == 3
	% 	h_dsm_np += delta1 * kron(CX3Dnp - CY3Dnp, kron(sigma_x, antidiag)) + delta2 * kron(CX3Dnp - CY3Dnp, kron(sigma_x, antidiag))
	end

	clear C*p S*p;
	[states_np,energy_eigenvalues_np] = eig(h_dsm_np);
	clear h_dsm_np;

	energy_eigenvalues_np = diag(energy_eigenvalues_np);
	no_of_orbitals = 2 * (2 * q + 1) * m * n * p

	no_of_zero_en_states = 0;
	for a = 1:no_of_orbitals
		if(abs(energy_eigenvalues_np(a)) < 10^-10)
			no_of_zero_en_states += 1;
		end
	end
	no_of_zero_en_states

	highest_state_np1 = states_np(:,no_of_orbitals/2);
	highest_state_np2 = states_np(:,(no_of_orbitals/2) - 1);
	% if q == 1
	% 	highest_state_np3 = states_np(:,(no_of_orbitals/2) - 2);
	% end
	highest_state_np = zeros(no_of_orbitals);

	% for a = 0:(level-1)
	% 	highest_state_np += abs(states_np(:,no_of_orbitals/2 -a)).^2;
	% end
	density_np = zeros(n,m,p);
	% In 3D, (i,k,l) is mapped to (l-1)*m*n + (k-1) * m + i, where l varies from 1 to p, k from 1 to n, i from 1 to m
	for a = 1:m
		for b = 1:n
			for c = 1:p
				for d = 0:(2* (2 * q + 1) - 1) %Number of internal states are 2 * (2*q + 1), where q is the spin (i.e. j)
					% a,b are intentionally swapped
					density_np(b,a,c) += abs(highest_state_np1(2*(2 * q + 1) *((c-1)*m*n + (b-1)*m + a) - d))^2 + abs(highest_state_np2(2*(2 * q + 1) *((c-1)*m*n + (b-1)*m + a) - d))^2;
					% density_np(b,a,c) += abs(highest_state_np(2*(2 * q + 1) *((c-1)*m*n + (b-1)*m + a) - d))^2;
				end
			end
		end
	end
	clear highest_state_np;
	max_density = max(max(max(density_np)))
	length(density_np)
	density_np = density_np/max_density;
	density_np(abs(density_np)<cutoff) = 0;

	[xx,yy,zz] = meshgrid(linspace(1,m,m),linspace(1,n,n),linspace(1,p,p));


	figure(1);
	scatter(linspace(1, no_of_orbitals, no_of_orbitals), energy_eigenvalues_np,'.');
	xlabel('$m$', 'fontsize', 20);
	ylabel('$E$', 'fontsize', 20);
	% title("Allowed energies for OBC", 'fontsize', 20)
	set(gca, "linewidth", 2, "fontsize", 20);
	axis tight;
	box on;

	title2 = strcat('OBC,','$E = ',num2str(energy_eigenvalues_np(no_of_orbitals/2)),'$')

	figure(2)
	scatter3(xx(:),yy(:),zz(:),100*density_np(:),density_np(:),"fill");
	colormap(rainbow);
	colorbar('fontsize',20);
	set(gca, "linewidth", 2, "fontsize", 20);
	title(title2, 'fontsize', 30);
	xlabel('$x$','fontsize', 20);
	ylabel('$y$','fontsize', 20);
	zlabel('$z$','fontsize', 20);
	xlim([1 m]);
	ylim([1 n]);
	zlim([1 p]);
	hold on
	plot3([1,m],[1,1],[1,1],'linewidth',2,'k')
	plot3([m,m],[n,n],[1,1],'linewidth',2,'k')
	plot3([m,1],[n,n],[1,1],'linewidth',2,'k')
	plot3([1,1],[n,1],[1,1],'linewidth',2,'k')
	plot3([1,m],[1,1],[p,p],'linewidth',2,'k')
	plot3([m,m],[1,n],[p,p],'linewidth',2,'k')
	plot3([m,1],[n,n],[p,p],'linewidth',2,'k')
	plot3([1,1],[n,1],[p,p],'linewidth',2,'k')
	plot3([1,1],[1,1],[1,p],'linewidth',2,'k')
	plot3([1,1],[n,n],[1,p],'linewidth',2,'k')
	plot3([m,m],[1,1],[1,p],'linewidth',2,'k')
	plot3([m,m],[n,n],[1,p],'linewidth',2,'k')

	toc;
	tic;
	cd saved_plots/strongtopo3D/energies
	mkdir(strcat('s',num2str(q)));
	cd(strcat('s',num2str(q)));
	name1 = strcat('order',num2str(topological_order),'OBC','lx',num2str(m),'ly',num2str(n),'lz',num2str(p));

	print(figure(1),'-dpdflatexstandalone',name1);

	system(strcat("pdflatex\t",name1));
	system('rm *.log *.aux');
	
	cd ../../../..

	cd saved_plots/strongtopo3D/density_plots
	mkdir(strcat('s',num2str(q)));
	cd(strcat('s',num2str(q)));

	print(figure(2),'-dpdflatexstandalone',name1);

	system(strcat("pdflatex\t",name1));
	system('rm *.log *.aux');
	
	cd ../../../..
	toc;
end 