classdef MagFEM_v_Adjoint < matlab.mixin.Copyable
    %   This class enable to create a class for for adjoint derivative from
    %   FEM
    properties
        fem                                        % the MagFEM object
        prmtrs=struct('name',{},'v_type',{})       % gradient parameters
        n_p=0                                        % number of gradient parameters
        dp
        dg
    end
    
    methods
        function obj = MagFEM_v_Adjoint(FEM_2D,n_p)
            % Construct an instance of this class
            obj.fem=FEM_2D;
            z_g=zeros(FEM_2D.n_n,n_p);z_m=zeros(FEM_2D.n_e,n_p);
            obj.dp={z_g,z_g,z_m,z_m,[]};
            assembling(obj);
        end
        
        function specifyGeometryVariation(obj,edg,vctr,name)
            % Specify the edge edg with the direction vector vctr and the name of
            % variation
            obj.n_p=obj.n_p+1;
            obj.prmtrs(obj.n_p).name=name;
            obj.prmtrs(obj.n_p).v_type='gmtr';
            p=unique(obj.fem.e(1:2,ismember(obj.fem.e(5,:),edg)));
            
            X=obj.fem.p(:,p)';
            v=vctr(X(:,1),X(:,2));
            obj.dp{1}(p,obj.n_p)=v{1};
            obj.dp{2}(p,obj.n_p)=v{2};
        end
        function specifyMaterialVariation(obj,rgn,fun,type,name)
            % Specify the region rgn, the type (J, nu) and the name of variation
            obj.n_p=obj.n_p+1;
            obj.prmtrs(obj.n_p).name=name;
            obj.prmtrs(obj.n_p).v_type='mtrl';
            p=ismember(obj.fem.t(end,:),rgn);
            if strcmp(type,'J')
                obj.dp{3}(p,obj.n_p)=fun;
            elseif strcmp(type,'nu')
                obj.dp{4}(p,obj.n_p)=1;%%%%% to do
            else
                error('Undifined property type')
            end
        end
        
        function[Jdx,Jdy,Jdj,Jdu]=currentDen(obj, x,ti)
            if nargin<3
                if isempty(x),    ti=1:obj.n_e;
                else,ti= pointLocation(triangulation(obj.fem.t(1:3,:)', obj.fem.p.'),x);
                end
            end
            ne=length(ti);
            Jdx=sparse(ne,obj.fem.n_n);
            Jdy=Jdx;
            Jdu=sparse(ne,obj.fem.n_u);
            Jdj=sparse(ne,obj.fem.n_e);
            idx = sub2ind(size(Jdj), 1:ne, ti');
            Jdj(idx)=1;
        end
       
        function assembling(obj)
            [gp,w]=gaussquad(obj.fem.NbQdPt);
            w=w';
            [lambda,gradx,grady]=obj.fem.base_func(gp);
            coordinates=obj.fem.p.';
            elements=obj.fem.t(1:end-1,:).';
            m=obj.fem.Mtrl;Js=obj.fem.j;Hr=obj.fem.hr_fun;   nuB=obj.fem.nl_fun;
            nF=size(gradx,2);nE=obj.fem.n_e;nN=obj.fem.n_n;unk=~obj.fem.bndr;
            if size(gradx,1)==1,nG=1;w1=.5;else,w1=w;nG=size(gradx,1);end
            % compute jac_1
            e=1e-100;z=zeros(1,2,1,6);
            z1=z;z1([1 8])=e*1j;
            z2=z;z2([3 10])=e*1j;
            z3=z;z3([5 12])=e*1j;
            c1 = coordinates(elements(:,1),:)+z1;
            d21 = coordinates(elements(:,2),:) - c1+z2;
            d31 = coordinates(elements(:,3),:) - c1+z3;
            
            dJac=d21(:,1,:,:).*d31(:,2,:,:)-d21(:,2,:,:).*d31(:,1,:,:);
            
            F=0;
            if  ~isempty(Hr)       % Magnet
                H=zeros(nE,2);
                for i=1:length(Hr)
                    H(m==200+i,:)=Hr{i};
                end
                F=F+w1*(gradx.*sum(d31.*H,2)-grady.*sum(d21.*H,2));
            end
            if ~isempty(Js)          % Copper  
                J=zeros(nE,1);
                for i=1:length(Js)
                    J(m==i)=Js{i};
                end
                F=F+w*lambda.*J.*dJac;
                djF=w*lambda.*real(dJac(:,1)).*obj.fem.mu0;
            end
            F=F.*obj.fem.mu0;
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
            
            U=obj.fem.u(elements);
            if ~isempty(nuB)          % Iron
                irnEl   = find(m>100 & m<200);
                xx=permute(gradx(:,id1).*gradx(:,id2),[3 2 1]);
                xy=permute(gradx(:,id1).*grady(:,id2),[3 2 1]);
                yx=permute(grady(:,id1).*gradx(:,id2),[3 2 1]);
                yy=permute(grady(:,id1).*grady(:,id2),[3 2 1]);
                
                vKn=(xx.*a(irnEl,:,:,:)+(xy+yx).*b(irnEl,:,:,:)+yy.*c(irnEl,:,:,:))./dJac(irnEl,:,:,:);
                Ku=permute(sum(reshape(vKn,[],nF,nF,nG,6).*U(irnEl,:),2),[1 3 4 5 2]);
                
                nB=sqrt(sum(Ku.*U(irnEl,:),2)./dJac(irnEl,:,:,:));
                nu=zeros(size(nB));
                for i=1:length(nuB)
                    id=m(irnEl)==100+i;
                    [nu(id,1,:),~]=nuB{i}(nB(id,1,:));
                end
                vKnn=sum(permute(w1,[1 3 2]).*nu.*vKn,3);
                
                vK(irnEl,:,:,:)=vKnn;
            end
            
            g=permute(imag(permute(sum(reshape(vK,[],nF,nF,6).*U,2),[1 3 2 4])-F)/e,[2 4 1 3]);
            
            iK=elements(:,repmat(1:nF,1,3))';
            jK=elements(:,reshape(repmat(1:3,nF,1),1,[]))';
            gx=reshape(g(:,1:3,:),[],1);
            gy=reshape(g(:,4:6,:),[],1);
            dxg=sparse(iK,jK,gx,nN,nN);
            dyg=sparse(iK,jK,gy,nN,nN);
            djg=sparse(elements,repmat(1:nE,1,nF),-djF,nN,nE);
            dxg=dxg(unk,:);dyg=dyg(unk,:);
            djg=djg(unk,:);
            [vK,iK,jK,b,NL_data]=obj.fem.initialization();
            if ~isempty(nuB)
                [~,dug]=obj.fem.setting_problem(obj.fem.u(unk),vK,iK,jK,b,NL_data);
            else
                dug=sparse(iK,jK,vK,nN,nN);
                dug=dug(unk,unk);
            end
            obj.dg={dxg,dyg,djg,dug};
        end
        
        function [Bx,By,dxBx,dyBx,dxBy,dyBy,duBx,duBy,dpBx,dpBy]=fluxDen(obj, x)
            coordinates=obj.fem.p.';
            ti= pointLocation(triangulation(obj.fem.t(1:3,:)', real(coordinates)),real(x));
            nE=length(ti);nN=obj.fem.n_n;
            unk=~obj.fem.bndr;
            elements=obj.fem.t(1:end-1,ti).';
            % compute jac_1
            e=1e-100;z=zeros(1,2,1,6);
            z1=z;z1([1 8])=e*1j;
            z2=z;z2([3 10])=e*1j;
            z3=z;z3([5 12])=e*1j;
            c1 = coordinates(elements(:,1),:)+z1;
            d21 = coordinates(elements(:,2),:) - c1+z2;
            d31 = coordinates(elements(:,3),:) - c1+z3;
            dJac=d21(:,1,:,:).*d31(:,2,:,:)-d21(:,2,:,:).*d31(:,1,:,:);
            
            X=x-c1;
            xr=[(d31(:,2).*X(:,1)-d31(:,1).*X(:,2)) ...
                (-d21(:,2).*X(:,1)+d21(:,1).*X(:,2))]./dJac;
            
            [~,gx,gy]=obj.fem.base_func(xr);
            
            Gradx=( d31(:,2).*gx-d21(:,2).*gy)./dJac;
            Grady=(-d31(:,1).*gx+d21(:,1).*gy)./dJac;
            
            U=obj.fem.u(elements);
            if length(ti)==1, U=U.';   end
            dxyBx=reshape( sum(imag(Grady)/e.*U,2),nE,3,[]);  
            dxyBy=reshape(-sum(imag(Gradx)/e.*U,2),nE,3,[]);
            
            iG=repmat(1:nE,1,3);  jG=elements(:,1:3).';
            iGu=repmat(1:nE,1,size(gx,2));jGu=elements.';
            dxBx=sparse(iG,jG,dxyBx(:,:,1),nE,nN);
            dxBy=sparse(iG,jG,dxyBy(:,:,1),nE,nN);
            dyBx=sparse(iG,jG,dxyBx(:,:,2),nE,nN);
            dyBy=sparse(iG,jG,dxyBy(:,:,2),nE,nN);
            duBx=sparse(iGu,jGu,real(Grady(:,:,1)),nE,nN);
            duBy=sparse(iGu,jGu,-real(Gradx(:,:,1)),nE,nN);
            duBx=duBx(:,unk);duBy=duBy(:,unk);
            
            c1=real(c1(:,:,1));
            d21=real(d21(:,:,1));
            d31=real(d31(:,:,1));
            d21=real(d21(:,:,1));
            dJac=real(dJac(:,1));
            X=x-c1+[e*1j 0];
            xr=[(d31(:,2).*X(:,1)-d31(:,1).*X(:,2)) ...
                (-d21(:,2).*X(:,1)+d21(:,1).*X(:,2))]./dJac;
            [~,gx1,gy1]=obj.fem.base_func(xr);
            X=x-c1+[0 e*1j];
            xr=[(d31(:,2).*X(:,1)-d31(:,1).*X(:,2)) ...
                (-d21(:,2).*X(:,1)+d21(:,1).*X(:,2))]./dJac;
            [~,gx2,gy2]=obj.fem.base_func(xr);
            gx=cat(3,gx1,gx2);
            gy=cat(3,gy1,gy2);
            Gradx=( d31(:,2).*gx-d21(:,2).*gy)./dJac;
            Grady=(-d31(:,1).*gx+d21(:,1).*gy)./dJac;
            dpBx=reshape(sum(imag(Grady)/e.*U,2),nE,[]);
            dpBy=reshape(-sum(imag(Gradx)/e.*U,2),nE,[]);
            
            Bx=reshape( sum(real(Grady(:,:,1)).*U,2),nE,[]);  
            By=reshape(-sum(real(Gradx(:,:,1)).*U,2),nE,[]);
        end
        
        function [E,dxE,dyE,djE,duE]=magneticEnergy(obj)
            [gp,w]=gaussquad(obj.fem.NbQdPt);
            w=w';
            lambda=obj.fem.base_func(gp);
            coordinates=obj.fem.p.';
            elements=obj.fem.t(1:end-1,:).';
            m=obj.fem.Mtrl;Js=obj.fem.j;Hr=obj.fem.hr_fun; 
            nE=obj.fem.n_e;nN=obj.fem.n_n;nF=size(lambda,2);
            % compute jac_1
            e=1e-100;z=zeros(1,2,1,6);
            z1=z;z1([1 8])=e*1j;
            z2=z;z2([3 10])=e*1j;
            z3=z;z3([5 12])=e*1j;
            c1 = coordinates(elements(:,1),:)+z1;
            d21 = coordinates(elements(:,2),:) - c1+z2;
            d31 = coordinates(elements(:,3),:) - c1+z3;
            
            dJac=d21(:,1,:,:).*d31(:,2,:,:)-d21(:,2,:,:).*d31(:,1,:,:);
            
            if  ~isempty(Hr)       % Magnet
                %%%%%%%%%%%%
            end
            if ~isempty(Js)          % Copper  
                J=zeros(nE,1);
                for i=1:length(Js)
                    J(m==i)=Js{i};
                end
                F=w*lambda.*J.*dJac ;
                djF=w*lambda.*real(dJac(:,1));
            end
            dxyF=permute(imag(F)/e,[2 4 1 3]);
            
            iK=elements(:,repmat(1:nF,1,3))';
            jK=elements(:,reshape(repmat(1:3,nF,1),1,[]))';
            Fx=reshape(dxyF(:,1:3,:),[],1);
            Fy=reshape(dxyF(:,4:6,:),[],1);
            
            u=.5*obj.fem.u.'; 
            % E=A'*f;
            dxE=u*sparse(iK,jK,Fx,nN,nN);
            dyE=u*sparse(iK,jK,Fy,nN,nN);
            djE=u*sparse(elements,repmat(1:nE,1,nF),djF,nN,nE);
            Eu=accumarray(elements(:),reshape(real(F(:,:,1)),[],1))';
            duE=.5*Eu(~obj.fem.bndr);
            lam=solving(obj,duE);
            if nargout==1
                E=   (dxE+lam*obj.dg{1})*obj.dp{1}...
                    +(dyE+lam*obj.dg{2})*obj.dp{2}...
                    +(djE+lam*obj.dg{3})*obj.dp{3};
            else
                E=.5*Eu*obj.fem.u;
            end
        end
        
        function lam=solving(obj,dydu)
            lam=-dydu/obj.dg{4};
        end
        
    end
end