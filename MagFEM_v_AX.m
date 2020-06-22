classdef MagFEM_v_AX < MagFEM_v
    %Reda_FEM_2D
    %   This class enable to create a class for axisymmetric
    
    methods
        function obj = MagFEM_v_AX(gmtr)
           % Construct an instance of this class
            obj@MagFEM_v(gmtr)
            obj.NbQdPt =6;
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
%                 H=zeros(nE,2);
%                 for i=1:length(Hr)
%                     H(m==200+i,:)=Hr{i};
%                 end
%                 F=F+accumarray(elements(:), reshape(w1*(gradx.*sum(d31.*H,2)-grady.*sum(d21.*H,2)) ,[],1));
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
            
            xx=permute(gradx(:,id1).*gradx(:,id2),[3 2 1]);
            xy=permute(gradx(:,id1).*grady(:,id2),[3 2 1]);
            yx=permute(grady(:,id1).*gradx(:,id2),[3 2 1]);
            yy=permute(grady(:,id1).*grady(:,id2),[3 2 1]);
            r_1=permute(w./(c1(:,1)+gp(:,1)'.*d21(:,1)+gp(:,2)'.*d31(:,1)),[1 3 2]);
            
            % matrices
            vK=reshape((sum(r_1.*(xx.*a+(xy+yx).*b+yy.*c),3)./dJac)',[],1);
                        
            % compute indices of sparse matrix
            iK=reshape(elements(:,id1)',[],1);
            jK=reshape(elements(:,id2)',[],1);
            
            
            if ~isempty(nuB)          % Iron
                error
            else,            NL_data=[];
            end
       end
       
        function [BB,AA,xx,yy,m,M]=computeDist(obj,N)
            x=gaussquad(N);
            [~,gradx,grady]=obj.base_func(x);
            %if size(gradx,1)==1, N=1;x=mean(x);end
            
            coordinates=obj.p.';
            elements=obj.t(1:end-1,:).';
            % compute jac_1
            c1 = coordinates(elements(:,1),:);
            d21 = coordinates(elements(:,2),:) - c1;
            d31 = coordinates(elements(:,3),:) - c1;
            dJac=d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
            
            x=permute(x,[3 2 1]);
            X=c1+x(1,1,:).*d21+x(1,2,:).*d31;
            r_1=1./X(:,1,:);
            X=reshape(permute(X,[2 1 3]),2,[]).';
            gx=permute(gradx,[3 2 1]);
            gy=permute(grady,[3 2 1]);
            Gradx=r_1.*( d31(:,2).*gx-d21(:,2).*gy)./dJac;
            Grady=r_1.*(-d31(:,1).*gx+d21(:,1).*gy)./dJac;
            
            U=obj.u(elements);
            Bx=sum(Gradx.*U,2);
            By=sum(Grady.*U,2);
            nB=sqrt(Bx.^2+By.^2);
            
            m=min(obj.p,[],2);M=max(obj.p,[],2);
            dx=sqrt(prod(M-m)/(2*size(obj.t,2)))/(2*N);
            [xx,yy] = meshgrid(m(1):dx:M(1),m(2):dx:M(2));
            F = scatteredInterpolant(X,nB(:),'natural','nearest');         BB=F(obj.p');
            AA= griddata(obj.p(1,:),obj.p(2,:),obj.u./obj.p(1,:).',xx,yy,'natural');AA(isnan(AA))=0;
            
        end
        
         function [Br,Bz,ti]=fluxDen(obj,x)
            
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
                Br=-Grady*U./x(:,1);  Bz=Gradx*U./x(:,1);
            else
                Br=-sum(Grady.*U,2)./x(:,1);  Bz=sum(Gradx.*U,2)./x(:,1);
            end
         end
        
       
    end
end

