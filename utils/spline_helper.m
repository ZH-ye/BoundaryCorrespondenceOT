classdef spline_helper < handle

    methods (Access = public)
        function this= spline_helper(use)
            this.order = [3,3]; % cubic % different order along each direction is not supported
            this.size_ctrl = [25,25]; % 25 control points
            this.uv_size=[20,20];% 2* 20 curves
            this = setup(this);
            this.u_handles=[];
            this.v_handles=[];
            this.surf_handle=[];
            this.colorbar_handle=[];

            this.patch_u_handle=[];

            this.patch_v_handle=[];
                        global use_export_setting;
            if isempty(use_export_setting) || ~ use_export_setting
            this.style_u = {'Color','red','PickableParts','none'};
            this.style_v = {'Color','blue','PickableParts','none'};
            else
            this.style_u = {'Color','red','PickableParts','none','LineWidth',1.25};
            this.style_v = {'Color','blue','PickableParts','none','LineWidth',1.25};
            end
            this.style_u_alt = {'Color',[1,0.5,0.5],'PickableParts','none'};
            this.style_v_alt = {'Color',[0.5,0.5,1],'PickableParts','none'};    
            this.gamma = 0.0001;
            this.nu_w = 0.0001;
            oldpath=path;
            path_this=mfilename('fullpath');
            sp=split(path_this,filesep);
            sp=sp(1:end-1);
            folder_this = strjoin(sp,filesep);
%             path(oldpath,[folder_this,filesep,'..',filesep,'BFIM']);
%             addpath([folder_this,filesep,'rulebased']);
            if nargin<1
                use = true;
                this.valid = use;
            end
        end
        function delete(sh)
            %sh.clear_draw();
        end
        function ready=is_ready(sh)
            ready=all(sh.current_fp>0);
        end
        function set_points(sh,points)
            % renew id and save points data
            if all(size(points)==[1024,2])
                points = points';
            end
            sh.points = points;
            sh.current_fp = zeros(4,1,'int32');
            sh.id = tempname;
            sh.clear_draw();

        end
        function save_coeff(sh,filename)
            mx = sh.size_ctrl(1);
            my = sh.size_ctrl(2);
            mat_coeff = zeros(mx*2,my);
            mat_coeff(1:mx,:)=sh.coeff(1,:,:);
            mat_coeff(mx+1:end,:)=sh.coeff(2,:,:);
            save(filename,'mat_coeff','-ascii');

        end
        function set_coeff(sh,coeff)
            sh.coeff = coeff;
            sh.sp = spmak({sh.knotx,sh.knoty},sh.coeff);
            if ~isempty(sh.on_coeff_change_callback)
                sh.on_coeff_change_callback();
            end
        end
        function load_coeff(sh,filename,tform)
            m_ptrfCoff = load(filename);
            [M,N] = size(m_ptrfCoff);
            u_basis_number = M/2;
            v_basis_number = N;
            sh.size_ctrl=[u_basis_number,v_basis_number];
            sh.setup();
            sh.coeff = zeros(2,u_basis_number,v_basis_number);
            sh.coeff(1,:,:)= m_ptrfCoff(1:u_basis_number,:);
            sh.coeff(2,:,:)= m_ptrfCoff(u_basis_number+1:end,:);
            if nargin==3
                [x,y] = tform.transformPointsInverse(sh.coeff(1,:,:),sh.coeff(2,:,:));
            sh.coeff(1,:,:) =x;
            sh.coeff(2,:,:) =y;
            end
            sh.sp = spmak({sh.knotx,sh.knoty},sh.coeff);
            if ~isempty(sh.on_coeff_change_callback)
                sh.on_coeff_change_callback();
            end
        end
        function reset_drawing(sh,axis,fp)
            init_map(sh,fp);
            sh.do_drawing_uv(axis,'uv');
        end
        function update_drawing(sh,axes,method,color)
            if ~all(sh.current_fp>0)
%                 sh.switch_alt(fp);
                return;
            end
%             Name = getfilename(sh,sh.current_fp);
%             if exist(Name, 'file') ~= 2
%                 return;
%             end
            if method(1)=='u' || method(1)=='U'
                sh.do_drawing_uv(axes,color);
            else
                if method(1)=='h' || method(1)=='H'
                do_drawing_heatmap(sh,axes,color)
                else
                    sh.clear_draw();
            end
            end

        end
        function try_update_spline(sh,fp,iter)
            if ~all(fp>0)
                return;
            end
            sh.current_fp = sort(reshape(fp,4,1));
            Name = getfilename(sh,sh.current_fp);
            if exist(Name, 'file') == 2
                load_coeff(sh,Name);
                sh.use_alt=false;
            else
                init_map(sh,fp);
                refine(sh,iter);
            end
        end
        function agree = check_coeff(sh,fp,coeff)
            if ~all(fp>0)
                agree = false;
                return;
            end
                     % check if fp and coeff matches
            coeff_fp = coeff(:,[1,end],[1,end]);
            coeff_fp = reshape(coeff_fp,2,[]);
            valid = true;
            for index=1:4
                i = fp(index);
                if valid
                    p = sh.points(:,i);
                    dist = min(sum(abs(coeff_fp-p),1));
                    if dist>=1e-3
                        valid = false;
                    end
                end
            end
            agree = valid;

        end
        function try_draw_coeff(sh,axis,fp,coeff,iter,form,color)
                      % try to load temp files
            if ~all(fp>0)
                sh.switch_alt(fp);
                return;
            end
            t = coeff == 0;
            if all(t(:))
                try_draw(sh,axis,fp,iter,form,color);
                return;
            end
%             if ~ sh.check_coeff(fp,coeff)
%                 try_draw(sh,axis,fp,iter,form,color);
%                 return;
%             end
%    
            
            sh.current_fp = sort(reshape(fp,4,1));

            sh.set_coeff(coeff);
            
            % drawing
            update_drawing(sh,axis,form,color);
        end
        function try_draw(sh,axis,fp,iter,form,color)

          % try to load temp files
            if ~all(fp>0)
                sh.switch_alt(fp);
                return;
            end
            try_update_spline(sh,fp,iter);
 
            % drawing
            update_drawing(sh,axis,form,color);

        end
        function try_draw_uv(sh,axis,fp,iter,color)
            if nargin ==4
                color = 'default';
            end
            try_draw(sh,axis,fp,iter,'uv',color)
        end
        function draw_heatmap(sh,axis,fp,iter,color)
        end
        function init_map(sh,fp,coeff)
            if ~all(fp>0)
                return;
            end
            sh.current_fp = sort(reshape(fp,4,1));
            fp = sh.current_fp;
            if nargin == 2
                fit_four_boundary(sh);
            else
                if nargin == 3 && ~ all(coeff ==0)
                    sh.set_coeff(coeff);
                end
            end
        end
        function refine(sh,iter)
            if ~all(sh.current_fp>0)
                return;
            end
            call_parametrization(sh,iter);
        end
        function refine_uv(sh,axis,iter)
            if ~all(sh.current_fp>0)
                sh.switch_alt(fp);
                return;
            end
            refine(sh,iter);

            sh.use_alt = false;
            sh.do_drawing(axis);
        end
        function switch_alt(sh,fp)
            if all(fp == sh.current_fp) && all(fp>0)
                should_use_alt = false;
            else
                should_use_alt = true;
            end
            if sh.use_alt ~= should_use_alt
                sh.use_alt = should_use_alt;
                sh.do_switch_style();
            end
        end
        function fit_four_boundary(sh)
            mx = sh.size_ctrl(1);
            my = sh.size_ctrl(2);
            sh.set_coeff(fit_boundary(mx,my,sh.knotx,sh.knoty,sh.points,sh.current_fp));
        end
        function call_parametrization(sh,iter)
            % call external program to update coefficients
            if iter ==0
                return
            end
            NAME = sh.getfilename(sh.current_fp);
            cornerpoints= sh.points(:,sort(sh.current_fp))';
            cc = cornerpoints;
            fixed = [-1,-1;1,-1;1,1;-1,1];
            dist_Min= 10000;
            for i = 1:4
                c = [cornerpoints(i+1:end,:);cornerpoints(1:i,:)];
                dd = fixed-c;
                dist = sum(dd(:).^2);
                if dist <dist_Min
                    dist_Min = dist;
                    cc = c;
                end
            end
            cornerpoints=cc;
            curv_length = to_curvature_length(cc');
            outerangle = curv_length(1,:);
            global use_trans;
            if ~isempty(use_trans)
                if use_trans
                    if any(outerangle<pi/3)
                        tform = fitgeotrans(fixed,fixed,'affine');
                    else
                        tform = fitgeotrans(cornerpoints,fixed,'affine');
                    end
                else
                    tform = fitgeotrans(fixed,fixed,'affine');
                    
                end
            else
                if any(outerangle<pi/3)
                    tform = fitgeotrans(fixed,fixed,'affine');
                else
                    tform = fitgeotrans(cornerpoints,fixed,'affine');
                end
            end

 
            
            %tform = fitgeotrans(cornerpoints,fixed,'affine');

            [x,y] = tform.transformPointsForward(sh.coeff(1,:,:),sh.coeff(2,:,:));
            sh.coeff(1,:,:) =x;
            sh.coeff(2,:,:) =y;
            sh.save_coeff(NAME);


            path_this=mfilename('fullpath');
            sp=split(path_this,filesep);
            sp=sp(1:end-1);
            folder_this = strjoin(sp,filesep);
            chars = @(n) convertStringsToChars(string(n));
            mx = sh.size_ctrl(1);
            my = sh.size_ctrl(2);
            order = sh.order(1);
            opt0 = [chars(order+1),' ',chars(mx-order-1),' ',chars(my-order-1),...
                ' 1 0.0 1.0 '];
            opt1 = [chars(sh.gamma),' ',chars(sh.nu_w)];
%             command = [fullfile(folder_this,'ImplicitParametrization.exe'),...
%                 ' ',opt0,opt1,' ',chars(iter),' ',NAME,' ',NAME];
            ld_library_path = strcat(getenv('LD_LIBRARY_PATH'),':',matlabroot,'/extern/bin/glnxa64');

            command = [fullfile(folder_this,'ImplicitParametrization'),...
                ' ',opt0,opt1,' ',chars(iter),' ',NAME,' ',NAME,' ',sh.session_name];
            %fprintf('%s\n',command);
            %systemcpplibpath = '/usr/lib/x86_64-linux-gnu/libstdc++.so.6';
            % the libstdc++ shipped with matlab is older than the system one
%             system(['export LD_PRELOAD=',systemcpplibpath,' ; ',...
%                 'export LD_LIBRARY_PATH=',ld_library_path,' ; ',command]);
            %system([ 'export LD_LIBRARY_PATH=',ld_library_path,' ; ',command]);
            mexParametrization(order+1,mx-order-1,my-order-1,1,0.0,1.0,sh.gamma,sh.nu_w,iter,NAME,NAME);
            sh.load_coeff(NAME,tform);
            delete(NAME);

        end
    end
    methods (Access = public)
        function name = getfilename(sh,fp)
            % return a file name based on sh.id and fp
            chars = @(n) convertStringsToChars(string(n));
            fp = sort(fp);
            name = [sh.id,sh.session_name,'_',chars(fp(1)),...
                '_',chars(fp(2)),...
                '_',chars(fp(3)),...
                '_',chars(fp(4))...
                ]; 
        end
        function clear_temp_file(sh,id)
            % clear temp file started with id (todo)
        end

        function sh = setup(sh)
            % set up knots data and staff
            mx = sh.size_ctrl(1);
            my = sh.size_ctrl(2);
            
            x0 = zeros(1,sh.order(1));
            x1 = ones(1,sh.order(1));
            y0 = zeros(1,sh.order(2));
            y1 = ones(1,sh.order(2));
            sh.knotx = [x0,linspace(0,1,mx-sh.order(1)+1),x1];
            sh.knoty = [y0,linspace(0,1,my-sh.order(2)+1),y1];
%           order = 3,3:
%             sh.knotsx =[0,0,0,linspace(0,1,mx-2),1,1,1];
%             sh.knotsy =[0,0,0,linspace(0,1,my-2),1,1,1];
        end
        function status = get_status(sh,strict)
            if nargin <2
                strict = false;
            end
            num_p = 200;
%             num_p = 1000;
            [pX,pY]= meshgrid(linspace(0,1,num_p),linspace(0,1,num_p));
            corner_size = 0.001;
            corners =...
            (mod(pX-1+0.5,1)-0.5).^2+(mod(pY-1+0.5,1)-0.5).^2<corner_size^2;
            mu_data=reshape(abs(spline_distortion_mu(sh.sp,[(pX(:))';(pY(:))'])),size(pX));
            mu_data_no_corners = mu_data(~corners);
            [jac_data,find_flip]=spline_distortion_Jac(sh.sp,[(pX(:))';(pY(:))']);
            jac_data=reshape(jac_data,size(pX));
            jac_data_no_corners = jac_data(~corners);
            if strict && ~find_flip && max(mu_data_no_corners(:))>0.90
                num_p = 200;
                for ii = 1:10
                    startx = 0.1*(ii-1);
                    for jj = 1:10
                        starty = 0.1*(jj-1);
                        [pX,pY]= meshgrid(linspace(startx,startx+0.1,num_p),linspace(starty,starty+0.1,num_p));
                        [~,find_flip]=spline_distortion_Jac(sh.sp,[(pX(:))';(pY(:))']);
                        if find_flip
                            fprintf('flipping!!\n');
                            break;
                        end
                    end
                    
                end
                
            end

            status.find_flip = find_flip;
            mean_jac = mean(jac_data(:));
            log_jac_data=log(abs(jac_data/mean_jac));
            mean_jac_no_corners = mean(jac_data_no_corners(:));
            log_jac_data_no_corners=log(abs(jac_data_no_corners/mean_jac_no_corners));
            status.min_log_jac = min(log_jac_data(:));
            status.max_log_jac = max(log_jac_data(:));
            status.min_log_jac_no_corners = min(log_jac_data_no_corners (:));
            status.max_log_jac_no_corners = max(log_jac_data_no_corners(:));
            status.max_mu = max(mu_data(:));
            status.max_mu_no_corners = max(mu_data_no_corners(:));
            status.l2_jac = std(jac_data(:));
            status.avr_mu = mean(mu_data(:));
            
        end
        function [logJs,Jmax,Jmin] = get_log_Js(sh)
            [~,Jgrid] = spline_distortion_Jac_Int(sh.sp,[]);
            Jgrid = abs(Jgrid);
            Jmin = min(Jgrid(:));
            Jmax = max(Jgrid(:));
%             Jstd = std(Jgrid(:));
            Jmean = mean(Jgrid(:));
            Jmax = Jmax /Jmean;
            Jmin = Jmin /Jmean;
            logJs = log(Jmax/Jmin);
        end
        function do_drawing_heatmap(sh,axis,color)
                        colormap(axis,'jet');

            clear_draw(sh);
            num_p = 200;
            hold(axis,'on');
            [pX,pY]= meshgrid(linspace(0,1,num_p),linspace(0,1,num_p));
            data = fnval(sh.sp,[(pX(:))';(pY(:))']);
            xdata = reshape(data(1,:),size(pX));
            ydata = reshape(data(2,:),size(pY));
            switch color(1)
                case 'm'
                    % mu  
%                                 caxis([0,1]);

                    mu_data=reshape(abs(spline_distortion_mu(sh.sp,[(pX(:))';(pY(:))'])),size(pX));
                    mu_data=remove_corner(mu_data,pX,pY,0);

                    sh.surf_handle = surf(axis,xdata,ydata,zeros(size(xdata)),mu_data,...
                        'EdgeColor','none','FaceColor','interp');
                    uistack(sh.surf_handle,'bottom');
                    sh.colorbar_handle= colorbar(axis);
                case 'l'
                    % log jacobian
%                                                     caxis([-2,2]);

                    jac_data=reshape(spline_distortion_Jac_Int(sh.sp,[(pX(:))';(pY(:))']),size(pX));
                    mean_jac = mean(jac_data(:));
                    %neg = jac_data<1e-10;
                    jac_data=log(abs(jac_data/mean_jac));
%                     jac_data = jac_data - mean(jac_data(:));
% %                     jac_data(jac_data<-3)=-3;
% %                     jac_data(neg)=-3;
% %                     jac_data(jac_data>3)=3;
                    %mm = min(jac_data(:));
%                     if mm<-1
%                         mm = -1.5;
%                     end
                    %jac_data(jac_data<mm)=mm;
                    %jac_data(neg)=mm;
                    %jac_data(jac_data>3)=3;
                    sh.surf_handle = surf(axis,xdata,ydata,zeros(size(xdata)),jac_data,...
                        'EdgeColor','none','FaceColor','interp');
                    uistack(sh.surf_handle,'bottom');
                    sh.colorbar_handle= colorbar(axis);

                otherwise
                    % red and blue
                    %sh.u_handles=plot(axis,xdata',ydata',sh.style_u{:});
            end
        end
        function do_drawing(sh,axis)
            sh.do_drawing_uv(axis,'red_and_blue');
        end
        function do_drawing_uv(sh,axis,color)
            colormap(axis,'jet');
            % calculate uv curves and draw
%             delete(sh.sf_handle);
%             for uh = sh.u_handles
%                 delete(uh);
%             end
%             for vh = sh.v_handles
%                 delete(vh);
%             end
            if nargin<2
                ff=figure;
                axis = axes(ff);
                color = 'u';
            end
            clear_draw(sh);
            num_p = 200;
            hold(axis,'on');
%           u
            [pX,pY]= meshgrid(linspace(0,1,num_p),linspace(0,1,sh.uv_size(2)));
            %pX,pY: 200*20
            data = fnval(sh.sp,[(pX(:))';(pY(:))']);
            xdata = reshape(data(1,:),size(pX));
            ydata = reshape(data(2,:),size(pY));
            switch color(1)
                case 'm'
%                                                     caxis([0,1]);

                    % mu  
                    mu_data=reshape(abs(spline_distortion_mu(sh.sp,[(pX(:))';(pY(:))'])),size(pX));
                    mu_data=remove_corner(mu_data,pX,pY,0);
                    xdata=xdata';
                    ydata=ydata';
                    mu_data=mu_data';
                    pX=pX';
                    xdata = [xdata; zeros(1,size(pX,2))];
                    ydata = [ydata; Inf*ones(1,size(pX,2))];
                    mu_data = [mu_data; zeros(1,size(pX,2))];
                    sh.patch_u_handle = patch(axis,xdata,ydata,mu_data,...
                        'EdgeColor','interp','FaceColor','none');
                    uistack(sh.patch_u_handle,'bottom');
                    %sh.colorbar_handle= colorbar(axis);
                case 'l'
%                     caxis([-2,2]);

                    % log jacobian
                    jac_data=reshape(spline_distortion_Jac_Int(sh.sp,[(pX(:))';(pY(:))']),size(pX));
                    mean_jac = mean(jac_data(:));
                    neg = jac_data<1e-3;
                    jac_data=log(abs(jac_data/mean_jac));
%                     jac_data = jac_data - mean(jac_data(:));
%                     jac_data(jac_data<-3)=-3;
%                     jac_data(neg)=-3;
%                     jac_data(jac_data>3)=3;
                    mm = min(jac_data(:));
%                     if mm<-1
%                         mm = -1.5;
%                     end
                      
%                     jac_data(jac_data<mm)=mm;
%                     jac_data(neg)=mm;
                    xdata=xdata';
                    ydata=ydata';
                    jac_data=jac_data';
                    pX=pX';
                    xdata = [xdata; zeros(1,size(pX,2))];
                    ydata = [ydata; Inf*ones(1,size(pX,2))];
                    jac_data = [jac_data; zeros(1,size(pX,2))];
                    sh.patch_u_handle = patch(axis,xdata,ydata,jac_data,...
                        'EdgeColor','interp','FaceColor','none');
                    uistack(sh.patch_u_handle,'bottom');
                otherwise
                    % red and blue
                    sh.u_handles=plot(axis,xdata',ydata',sh.style_u{:});
            end
            %v
            [pX,pY]= meshgrid(linspace(0,1,sh.uv_size(1)),linspace(0,1,num_p));
            data = fnval(sh.sp,[(pX(:))';(pY(:))']);
            xdata = reshape(data(1,:),size(pX));
            ydata = reshape(data(2,:),size(pY));
            switch color(1)
                case 'm'
%                     caxis([0,1]);

                    % mu  
                    mu_data=reshape(abs(spline_distortion_mu(sh.sp,[(pX(:))';(pY(:))'])),size(pX));
                    mu_data=remove_corner(mu_data,pX,pY,0);
                    xdata = [xdata; zeros(1,size(pX,2))];
                    ydata = [ydata; Inf*ones(1,size(pX,2))];
                    mu_data = [mu_data; zeros(1,size(pX,2))];
                    %sh.surf_handle = surf(xdata,ydata,zeros(size(xdata)),mu_data,'EdgeColor','interp',...
                    %    'FaceColor','none');
                    sh.patch_v_handle = patch(axis,xdata,ydata,mu_data,...
                        'EdgeColor','interp','FaceColor','none');
                    uistack(sh.patch_v_handle,'bottom');
                    sh.colorbar_handle= colorbar(axis);

                case 'l'
%                     caxis([-2,2]);

                    % log jacobian
                    jac_data=reshape(spline_distortion_Jac_Int(sh.sp,[(pX(:))';(pY(:))']),size(pX));
                    mean_jac = mean(jac_data(:));
                    neg = jac_data<1e-3;
                    jac_data=log(abs(jac_data/mean_jac));
                    jac_data = jac_data - mean(jac_data(:));
%                     jac_data(jac_data<-3)=-3;
%                     jac_data(neg)=-3;
%                     jac_data(jac_data>3)=3;
                    mm = min(jac_data(:));
%                     if mm<-1
%                         mm = -1.5;
%                     end
                    %mm = min(jac_data(:));
                    jac_data(jac_data<mm)=mm;
                    jac_data(neg)=mm;
                    xdata = [xdata; zeros(1,size(pX,2))];
                    ydata = [ydata; Inf*ones(1,size(pX,2))];
                    jac_data = [jac_data; zeros(1,size(pX,2))];
                    sh.patch_v_handle = patch(axis,xdata,ydata,jac_data,...
                        'EdgeColor','interp','FaceColor','none');
                    uistack(sh.patch_v_handle,'bottom');
                    sh.colorbar_handle= colorbar(axis);
                                        
                otherwise
                    % red and blue
                    sh.v_handles=plot(axis,xdata,ydata,sh.style_v{:});

            end
%             sh.sf_handle= surface('XData',xdata,'YData',ydata,'ZData',zeros(size(xdata)),'FaceColor','none');

        end   
        function do_switch_style(sh)
            % switch style based on current state of use_alt 
            if sh.use_alt
                su = sh.style_u_alt;
                sv = sh.style_v_alt;
            else
                su = sh.style_u;
                sv = sh.style_v;
            end
            for uh = sh.u_handles
                set(uh,su{:});
            end
            for vh = sh.v_handles
                set(vh,sv{:});
            end
        end
        function clear_draw(sh)
            %todo
            for uh = sh.u_handles
                delete(uh);
            end
            for vh = sh.v_handles
                delete(vh);
            end
            delete(sh.patch_u_handle);
            sh.patch_u_handle=[];
            delete(sh.patch_v_handle);
            sh.patch_v_handle=[];
            delete(sh.surf_handle);
            sh.surf_handle=[];
            %delete(sh.colorbar_handle);
            if ~isempty(sh.colorbar_handle) && isvalid(sh.colorbar_handle)
                colorbar(sh.colorbar_handle,'off');
            end
            sh.u_handles=[];
            sh.v_handles=[];
        end

    end
    properties (Access = public)
        valid
        id
        order
        size_ctrl
        % input data
        points % 2*1024
        current_fp % 4*1 int32
        % spline data
        coeff % 2*mx*my
        knotx
        knoty
        sp % need curve fitting toolbox
        status
        % drawing properties
        uv_size
        u_handles
        v_handles
        surf_handle
        patch_u_handle
        patch_v_handle
        colorbar_handle
        use_alt
        
        style_uv
        style_heatmap
        
        style_u
        style_u_alt
        style_v
        style_v_alt
        % external properties
        para_exec
        gamma
        nu_w
        session_name
        on_coeff_change_callback=[]
    end
end
function data = remove_corner(data,pX,pY,defaultvalue,corner_size)
if nargin<5
    corner_size = 0.001;
end
corners =...
    (mod(pX-1+0.5,1)-0.5).^2+(mod(pY-1+0.5,1)-0.5).^2<corner_size^2;
data(corners)=defaultvalue;
end