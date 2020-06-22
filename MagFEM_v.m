classdef MagFEM_v < matlab.mixin.Copyable
    %Reda_FEM_2D
    %   This class enable to create a class for 2D palnar and axisymmetric
    %   FEM
    properties
        type                        % type of the analysis planar 'plan'/ axisymmetric 'axi'
        gmtr                        % geometry created using PDE toolbox
        p                           % mesh coordinate
        t                           % mesh connectivity
        e                           % edges infos
        n_e                         % number of elements
        n_n                         % number of nodes
        n_u                         % number of unknowns
        bndr                        % nodes on dirchlet's boundary
        j                           % current densities             (indices below 100)
        nl_fun                      % function of BH-curve          (indices below 200)
        hr_fun                      % The magentisation directions  (indices below 300)
        Mtrl                        % Indices of differnt material
        mu0=4*pi*1e-7
        NbQdPt =4                   % Number of points in gauss quadrature
        base_func=@ref_base_func_P2     % base functions
        u                           % solution of the FEM problem
    end
    
    methods
        function obj = MagFEM_v(gmtr)
           % Construct an instance of this class
            obj.gmtr=gmtr;
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            
        end
        function  showGeometry(obj)
            % Display the geometry
            h2=pdegplot(obj.gmtr,'FaceLabels','on')%,'EdgeLabels','on');
            set(h2,'Color','m', 'LineWidth',2)
            axis equal
            addgradient(gca, [0 0 1],[1 1 1]);  
            %set(pl,'FaceAlpha',.3)  % Make transparent
        end
        
        function meshing(obj,varargin)
            % Mesh the geometry
            model = createpde(1);
            geometryFromEdges(model,obj.gmtr);
            generateMesh(model,varargin{:});
            [obj.p,obj.e,obj.t]= meshToPet(model.Mesh);
            obj.n_n=size(obj.p,2);
            obj.n_e=size(obj.t,2);
            obj.Mtrl=zeros(obj.n_e,1);
            obj.u=zeros(obj.n_n,1);
        end
        function refineMesh(obj,varargin)
            [obj.p,obj.e,obj.t] = refinemesh(obj.gmtr,obj.p,obj.e,obj.t,varargin{:});
            obj.p = jigglemesh(obj.p,obj.e,obj.t,'Iter',100);
            obj.n_n=size(obj.p,2);
            obj.n_e=size(obj.t,2);
            obj.Mtrl=zeros(obj.n_e,1);
            obj.u=zeros(obj.n_n,1);
        end
        function  showMesh(obj,rgn)
            % Display the mesh
            if nargin>1
                el=ismember(obj.t(end,:),rgn);
                pdeplot(obj.p,obj.e,obj.t(:,el))
            else
                pdeplot(obj.p,obj.e,obj.t)
            end
            %set(h2,'Color','m', 'LineWidth',2)
            axis equal
            m=min(obj.p,[],2);M=max(obj.p,[],2);
            xlim([m(1) M(1)]+max(M-m).*[-.2 .2])
            ylim([m(2) M(2)]+max(M-m).*[-.2 .2])
            addgradient(gca, [0 0 1],[1 1 1]);
        end
        
        function applyBoundary(obj,bndr)
            % Apply dirichlet homogeneous boudary conditions
            indiceE=ismember(obj.e(5,:),bndr);
            temp=unique(obj.e(1:2,indiceE));
            obj.bndr=ismember(1:obj.n_n,temp).';
            obj.n_u=obj.n_n-length(temp);
        end
        function applyCurrentD(obj,rgn,J)
            % Apply current density J to region(s) rgn
            obj.j{end+1}=J;
            s=size(obj.j,2);
            for r=rgn
                obj.Mtrl(obj.t(end,:)==r)=s;
            end
        end
        function applyMagnet(obj,rgn,hr)
            % Apply coercitive field hr to region(s) rgn
            obj.hr_fun{end+1}=hr;
            s=size(obj.hr_fun,2)+200;
            for r=rgn
                obj.Mtrl(obj.t(end,:)==r)=s;
            end
        end
        function applyMaterial(obj,rgn,fun)
            % Apply ferromagnetic material to region(s) rgn
            obj.nl_fun{end+1}=fun;
            s=size(obj.nl_fun,2)+100;
            for r=rgn
                obj.Mtrl(obj.t(end,:)==r)=s;
            end
        end
        
        function set_basis_function(obj,P)
            if strcmp(P,'P1')
                obj.base_func=@ref_base_func_P1;
            end
        end
        
       function [vK,iK,jK,F,NL_data]=initialization(obj)
            
            [gp,w]=gaussquad(obj.NbQdPt);
            w=w';
            [lambda,gradx,grady]=obj.base_func(gp);
            coordinates=obj.p.';
            elements=obj.t(1:end-1,:).';
            m=obj.Mtrl;Js=obj.j;Hr=obj.hr_fun;   nuB=obj.nl_fun;
            nF=size(gradx,2);nE=obj.n_e;
            if size(gradx,1)==1,w1=.5;else,w1=w;end
            % compute jac_1
            c1 = coordinates(elements(:,1),:);
            d21 = coordinates(elements(:,2),:) - c1;
            d31 = coordinates(elements(:,3),:) - c1;
            
            dJac=d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
            
            F=0;
            if  ~isempty(Hr)       % Magnet
                H=zeros(nE,2);
                for i=1:length(Hr)
                    H(m==200+i,:)=Hr{i};
                end
                F=F+accumarray(elements(:), reshape(w1*(gradx.*sum(d31.*H,2)-grady.*sum(d21.*H,2)) ,[],1));
            end
            if ~isempty(Js)          % Copper  
                J=zeros(nE,1);
                for i=1:length(Js)
                    J(m==i)=Js{i};
                end
                F=F+accumarray(elements(:), reshape(w*lambda.*J.*dJac ,[],1));
            end
            F=F.*obj.mu0;
            a=sum(d31.*d31,2);
            b=-sum(d31.*d21,2);
            c=sum(d21.*d21,2);
            %compute gradient
            id1=repmat(1:nF,1,nF);
            id2=reshape(repmat(1:nF,nF,1),1,[]);
            
            xx=w1*(gradx(:,id1).*gradx(:,id2));
            xy=w1*(gradx(:,id1).*grady(:,id2));
            yx=w1*(grady(:,id1).*gradx(:,id2));
            yy=w1*(grady(:,id1).*grady(:,id2));
            
            % matrices
            vK=((xx.*a+(xy+yx).*b+yy.*c)./dJac);
                        
            % compute indices of sparse matrix
            iK=reshape(elements(:,id1)',[],1);
            jK=reshape(elements(:,id2)',[],1);
            
            
            if ~isempty(nuB)          % Iron
                irnEl   = find(m>100 & m<200);
                nu=ones(nE,1);
                [nu(irnEl),~]=nuB{1}(0);
                vK=vK.*nu;
                
                xx=permute(gradx(:,id1).*gradx(:,id2),[3 2 1]);
                xy=permute(gradx(:,id1).*grady(:,id2),[3 2 1]);
                yx=permute(grady(:,id1).*gradx(:,id2),[3 2 1]);
                yy=permute(grady(:,id1).*grady(:,id2),[3 2 1]);
                
                vKt=(xx.*a(irnEl)+(xy+yx).*b(irnEl)+yy.*c(irnEl))./dJac(irnEl);
                
                NL_data.K= vKt;
                NL_data.id=m(irnEl);
                NL_data.idsp=irnEl;
                NL_data.elements=elements(irnEl,:);
                NL_data.dJac=dJac(irnEl);
                NL_data.w=permute(w1,[1 3 2]);
            else,            NL_data=[];
                vK=vK.';
            end
       end
       
       function [R,dR]=setting_problem(obj,u,vK,iK,jK,b,NL_data)
           unk=~obj.bndr;nn=obj.n_n;obj.u(unk)=u;
           U=obj.u(NL_data.elements);
           vKn=NL_data.K;
           w=NL_data.w;
           fun=obj.nl_fun;
           
           nF=sqrt(size(vKn,2));nG=size(vKn,3);
           id1=repmat(1:nF,1,nF);
           id2=reshape(repmat(1:nF,nF,1),1,[]);
           Ku=permute(sum(reshape(vKn,[],nF,nF,nG).*U,2),[1 3 4 2]);
           KuTKu=Ku(:,id1,:).*Ku(:,id2,:)./NL_data.dJac;
           
           nB=sqrt(sum(Ku.*U,2)./NL_data.dJac);
           nu=zeros(size(nB));dnu=zeros(size(nB));
           for i=1:length(fun)
               id=NL_data.id==100+i;
               [nu(id,1,:),dnu(id,1,:)]=fun{i}(nB(id,1,:));
           end
           M=dnu./nB; M(isinf(M))=0;
           vKt=sum(w.*nu.*vKn,3);
           vJKt=sum(w.*M.*KuTKu,3);
           
           vK(NL_data.idsp,:)=vKt;
           vJK=vK;
           vJK(NL_data.idsp,:)=vKt+vJKt;
           R=sparse(iK,jK,vK.',nn,nn)*obj.u-b;
           R=R(unk);
           JK=sparse(iK,jK,vJK.',nn,nn);
           dR=JK(unk,unk);
       end
        
        function solving(obj,varargin)
            nn=obj.n_n;unk=~obj.bndr;
            [vK,iK,jK,b,NL_data]=obj.initialization();
            
            K=sparse(iK,jK,vK.',nn,nn);
            A=K(unk,unk)\b(unk);
            
            if ~isempty(NL_data)
                options = optimoptions('fsolve','SpecifyObjectiveGradient',true,...'CheckGradients',true,...
                    varargin{:});
                objFun=@(u)obj.setting_problem(u,vK,iK,jK,b,NL_data);
                
                A=fsolve(objFun,A,options);
            end
            obj.u(~obj.bndr)=A;
            
        end
        
        function E=magneticEnergy(obj) 
            [gp,w]=gaussquad(obj.NbQdPt);
            w=w';
            [lambda,gradx,grady]=obj.base_func(gp);
            coordinates=obj.p.';
            elements=obj.t(1:end-1,:).';
            m=obj.Mtrl;Js=obj.j;Hr=obj.hr_fun; 
            nE=obj.n_e;
            if size(gradx,1)==1,w1=.5;else,w1=w;end
            % compute jac_1
            c1 = coordinates(elements(:,1),:);
            d21 = coordinates(elements(:,2),:) - c1;
            d31 = coordinates(elements(:,3),:) - c1;
            
            dJac=d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
            
            F=0;
            if  ~isempty(Hr)       % Magnet
                H=zeros(nE,2);
                for i=1:length(Hr)
                    H(m==200+i,:)=Hr{i};
                end
                F=F+accumarray(elements(:), reshape(w1*(gradx.*sum(d31.*H,2)-grady.*sum(d21.*H,2)) ,[],1));
            end
            if ~isempty(Js)          % Copper  
                J=zeros(nE,1);
                for i=1:length(Js)
                    J(m==i)=Js{i};
                end
                F=F+accumarray(elements(:), reshape(w*lambda.*J.*dJac ,[],1));
            end
            
            E=0.5*obj.u.'*F;
        end
        
        function showSolution(obj)
            N=obj.NbQdPt;
            [BB,AA,xx,yy,m,M]=obj.computeDist(N);
            
            clf('reset')
            axis equal
            xlim([m(1) M(1)]+max(M-m).*[-.2 .2])
            ylim([m(2) M(2)]+max(M-m).*[-.2 .2])
            xlabel('x')
            ylabel('y')
            title('Flux density')
            addgradient(gca, [0 0 1],[1 1 1]);  
            hold on
            %trisurf(elements,X(:,1),X(:,2),BB,'LineStyle','none') 
            %contourf(xx,yy,bb,20,'LineStyle','none');
            if size(obj.t,1)==4, idp=1:3;else, idp=[1 4 2 5 3 6];end
            patch('Faces',obj.t(idp,:)','Vertices',obj.p',...
                'FaceVertexCData',BB,'FaceColor','interp','edgecolor','none');
            colormap('jet'),caxis([0 max(BB(:))]),cbar=colorbar('south');
            contour(xx,yy,AA,linspace(min(obj.u),max(obj.u),20),'k');
            h2=pdegplot(obj.gmtr);
            set(h2,'Color','m', 'LineWidth',2)
            hold off
            ax = gca;
            axpos = ax.Position;
            cbP=cbar.Position ;
            cbar.Position = [cbP(1)+cbP(3)/4  cbP(2)+cbP(4)/2  cbP(3)/2  cbP(4)/2];
            ax.Position = axpos;
            set(ax,'FontSize',14)
            pause(.001)
        end
        function [BB,AA,xx,yy,m,M]=computeDist(obj,N)
            x=gaussquad(N);
            [~,gradx,grady]=obj.base_func(x);
            if size(gradx,1)==1, N=1;x=mean(x);end
            
            coordinates=obj.p.';
            elements=obj.t(1:end-1,:).';
            % compute jac_1
            c1 = coordinates(elements(:,1),:);
            d21 = coordinates(elements(:,2),:) - c1;
            d31 = coordinates(elements(:,3),:) - c1;
            dJac=d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
            
            x=permute(x,[3 2 1]);
            X=reshape(permute(c1+x(1,1,:).*d21+x(1,2,:).*d31,[2 1 3]),2,[])';
            
            gx=permute(gradx,[3 2 1]);
            gy=permute(grady,[3 2 1]);
            Gradx=( d31(:,2).*gx-d21(:,2).*gy)./dJac;
            Grady=(-d31(:,1).*gx+d21(:,1).*gy)./dJac;
            
            U=obj.u(elements);
            Bx=sum(Gradx.*U,2);
            By=sum(Grady.*U,2);
            nB=sqrt(Bx.^2+By.^2);
            
            m=min(obj.p,[],2);M=max(obj.p,[],2);
            dx=sqrt(prod(M-m)/(2*size(obj.t,2)))/(2*N);
            [xx,yy] = meshgrid(m(1):dx:M(1),m(2):dx:M(2));
            F = scatteredInterpolant(X,nB(:),'natural','nearest');         BB=F(obj.p');
            AA= griddata(obj.p(1,:),obj.p(2,:),obj.u,xx,yy,'natural');AA(isnan(AA))=0;
            
        end
        
         function [Bx,By,ti]=fluxDen(obj,x)
            
            ti= pointLocation(triangulation(obj.t(1:3,:)', real(obj.p)'),real(x));
            coordinates=obj.p.';
            elements=obj.t(1:end-1,ti).';
            % compute jac_1
            c1 = coordinates(elements(:,1),:);
            d21 = coordinates(elements(:,2),:) - c1;
            d31 = coordinates(elements(:,3),:) - c1;
            dJac=d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
            
            X=x-c1;
            xr=[(d31(:,2).*X(:,1)-d31(:,1).*X(:,2)) ...
               (-d21(:,2).*X(:,1)+d21(:,1).*X(:,2))]./dJac;
                
            [~,gx,gy]=obj.base_func(xr);
            
            Gradx=( d31(:,2).*gx-d21(:,2).*gy)./dJac;
            Grady=(-d31(:,1).*gx+d21(:,1).*gy)./dJac;
            
            U=obj.u(elements);
            if length(ti)==1
                Bx=Grady*U;  By=-Gradx*U;
            else
                Bx=sum(Grady.*U,2);  By=-sum(Gradx.*U,2);
            end
         end
        
         function f_nu=fitBHcurve(obj,B,H,type)
            switch type
                case 'spline'
                    nany=(B==0);
                    nur=H./B*obj.mu0;
                    B(nany)=[];nur(nany)=[];
                    pp=pchip(B,nur);
                    [breaks,coefs,l,k,d] = unmkpp(pp);
                    dpp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
                    f_nu=@(b)deal(ppval(pp,b), ppval(dpp,b));
                    %f_nu=[pp,dpp];
                case 'tanh'
                    f=@(p)H.*(pi*4e-7*( 1.0 + p(1)*(1-tanh(p(2)*B).^p(3))))-B;
                    [pm,fm]=fminsearch(@(x)norm(f(x)),[1000,3,200]);
                    %
                    %
                case 'marrocco'
                    %
                    %
            end
         end
    end
end

