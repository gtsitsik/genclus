function  execution_time_on_artificial_data()
% -------------- NODES -----------
disp('=====NODES COMPARISON STARTED=====')

size_all= {...
    {30,20,10},...
    {50,10},...
    {10,50}...
    };
noise_level = 0.01;
sparsity_level_all = [0.85];
node_mult=2.^[0:5];
for k = 1:numel(node_mult)
    params.graph_tree(k) = graph_tree_root;
    for i = 1:numel(size_all)
        params.graph_tree(k).Children(i).slices_num = 3;
        params.graph_tree(k).Children(i).noise_level = 0.01;
        params.graph_tree(k).Children(i).sparsity_level = 0.85;
        for j = 1:numel(size_all{i})
            params.graph_tree(k).Children(i).Children(j).type = 'clique';
            params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j}*node_mult(k);
        end
    end
    params.graph_tree(k).labels = [];
    [~,params.graph_tree(k)] = create_graph(params.graph_tree(k));
    %     params.graph_tree(k)=cur_graph_tree;
end
% 
prms.embedding_method.GenClus.A_const = '+';
prms.embedding_method.GenClus.B_const = '+';
prms.embedding_method.GenClus.rho = 0;
prms.embedding_method.ComClus.beta = linspace(0.01,0.9145,3);
prms.embedding_method.ComClus.rho = linspace(0,0.16,3);
prms.embedding_method.ComClus.thres_inner = [1e-6];
prms.embedding_method.Symmetric_Richcom.structure = "true";
prms.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
prms.embedding_method.CMNC.structure = "true";
prms.embedding_method.CMNC.delta = 1;  

params.thres = [1e-6];
params.max_iters = [1000];
params.R = 3;%*2.^[1:4];;
params.M = 3;%*2.^[1:3:15];
params.L_type_ind = [1 2];
params.sample = 1:5;


figure('Name','Number of Nodes')
sel_ind = 2;
x_dat=node_mult*60;

ax = axes;
hold on
ax.YScale='log';
ax.XScale='log';
ylabel('Duration (sec)')
xlabel('Number of nodes (I)')
legend
grid

params.embedding_method=[];
params.embedding_method.GenClus=prms.embedding_method.GenClus;
durations=speed_calculate(params);
create_plot('GenClus',x_dat,durations.GenClus,sel_ind)

params.embedding_method=[];
params.embedding_method.ComClus=prms.embedding_method.ComClus;
durations=speed_calculate(params);
create_plot('ComClus',x_dat,durations.ComClus,sel_ind)

params.embedding_method=[];
params.embedding_method.Symmetric_Richcom=prms.embedding_method.Symmetric_Richcom;
durations=speed_calculate(params);
create_plot('Symmetric Richcom',x_dat,durations.Symmetric_Richcom,sel_ind)

params.graph_tree=params.graph_tree(1:3);
params.embedding_method=[];
params.embedding_method.CMNC=prms.embedding_method.CMNC;
durations=speed_calculate(params);
create_plot('CMNC',x_dat(1:3),durations.CMNC,sel_ind)


disp('=====NODES COMPARISON ENDED =====')


%---------------- VIEWS-------------
disp('=====VIEWS COMPARISON STARTED=====')
clear all

size_all= {...
    {30,20,10},...
    {50,10},...
    {10,50}...
    };
noise_level = 0.01;
sparsity_level_all = [0.85];
node_mult=2.^[0:8];
for k = 1:numel(node_mult)
    params.graph_tree(k) = graph_tree_root;
    for i = 1:numel(size_all)
        params.graph_tree(k).Children(i).slices_num = 3*node_mult(k);
        params.graph_tree(k).Children(i).noise_level = 0.01;
        params.graph_tree(k).Children(i).sparsity_level = 0.85;
        for j = 1:numel(size_all{i})
            params.graph_tree(k).Children(i).Children(j).type = 'clique';
            params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j};
        end
    end
    params.graph_tree(k).labels = [];
    [~,params.graph_tree(k)] = create_graph(params.graph_tree(k));
    %     params.graph_tree(k)=cur_graph_tree;
end
% 
prms.embedding_method.GenClus.A_const = '+';
prms.embedding_method.GenClus.B_const = '+';
prms.embedding_method.GenClus.rho = 0;
prms.embedding_method.ComClus.beta = linspace(0.01,0.9145,3);
prms.embedding_method.ComClus.rho = linspace(0,0.16,3);
prms.embedding_method.ComClus.thres_inner = [1e-6];
prms.embedding_method.Symmetric_Richcom.structure = "true";
prms.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
prms.embedding_method.CMNC.structure = "true";
prms.embedding_method.CMNC.delta = 1;  

params.thres = [1e-6];
params.max_iters = [1000];
params.R = 3;%*2.^[1:4];;
params.M = 3;%*2.^[1:3:15];
params.L_type_ind = [1 2];
params.sample = 1:5;


figure('Name','Number of Views')
sel_ind = 2;
x_dat=node_mult*9;

ax = axes;
hold on
ax.YScale='log';
ax.XScale='log';
ylabel('Duration (sec)')
xlabel('Number of views (K)')
legend
grid

params.embedding_method=[];
params.embedding_method.GenClus=prms.embedding_method.GenClus;
durations=speed_calculate(params);
create_plot('GenClus',x_dat,durations.GenClus,sel_ind)

params.embedding_method=[];
params.embedding_method.ComClus=prms.embedding_method.ComClus;
durations=speed_calculate(params);
create_plot('ComClus',x_dat,durations.ComClus,sel_ind)

params.embedding_method=[];
params.embedding_method.Symmetric_Richcom=prms.embedding_method.Symmetric_Richcom;
params.graph_tree=params.graph_tree;
durations=speed_calculate(params);
create_plot('Symmetric Richcom',x_dat,durations.Symmetric_Richcom,sel_ind)
params.graph_tree=params.graph_tree(1:6);
params.embedding_method=[];
params.embedding_method.CMNC=prms.embedding_method.CMNC;
durations=speed_calculate(params);
create_plot('CMNC',x_dat(1:6),durations.CMNC,sel_ind)


disp('=====VIEWS COMPARISON ENDED =====')


%------------- R --------------
disp('=====RANK COMPARISON STARTED=====')
clear all

size_all= {...
    {30,20,10},...
    {50,10},...
    {10,50}...
    };
noise_level = 0.01;
sparsity_level_all = [0.85];
k=1;
params.graph_tree(k) = graph_tree_root;
for i = 1:numel(size_all)
    params.graph_tree(k).Children(i).slices_num = 3;
    params.graph_tree(k).Children(i).noise_level = 0.01;
    params.graph_tree(k).Children(i).sparsity_level = 0.85;
    for j = 1:numel(size_all{i})
        params.graph_tree(k).Children(i).Children(j).type = 'clique';
        params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j}*2^2;
    end
end
params.graph_tree(k).labels = [];
[~,params.graph_tree(k)] = create_graph(params.graph_tree(k));
%     params.graph_tree(k)=cur_graph_tree;
% 
prms.embedding_method.GenClus.A_const = '+';
prms.embedding_method.GenClus.B_const = '+';
prms.embedding_method.GenClus.rho = 0;
prms.embedding_method.ComClus.beta = linspace(0.01,0.9145,3);
prms.embedding_method.ComClus.rho = linspace(0,0.16,3);
prms.embedding_method.ComClus.thres_inner = [1e-6];
prms.embedding_method.Symmetric_Richcom.structure = "true";
prms.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
prms.embedding_method.CMNC.structure = "true";
prms.embedding_method.CMNC.delta = 1;  

params.thres = [1e-6];
params.max_iters = [1000];
params.R = 3*2.^[0:6];;
params.M = 3;%*2.^[1:3:15];
params.L_type_ind = [1 2];
params.sample = 1:5;


figure('Name','Rank (R)')
sel_ind = 3;
x_dat=params.R;

ax = axes;
hold on
ax.YScale='log';
ax.XScale='log';
ylabel('Duration (sec)')
xlabel('Rank (R)')
legend
grid

params.embedding_method=[];
params.embedding_method.GenClus=prms.embedding_method.GenClus;
durations=speed_calculate(params);
create_plot('GenClus',x_dat,durations.GenClus,sel_ind)

params.embedding_method=[];
params.embedding_method.ComClus=prms.embedding_method.ComClus;
durations=speed_calculate(params);
create_plot('ComClus',x_dat,durations.ComClus,sel_ind)

params.embedding_method=[];
params.embedding_method.Symmetric_Richcom=prms.embedding_method.Symmetric_Richcom;
durations=speed_calculate(params);
create_plot('Symmetric Richcom',x_dat,durations.Symmetric_Richcom,sel_ind)

params.R=params.R(1:2);
params.embedding_method=[];
params.embedding_method.CMNC=prms.embedding_method.CMNC;
durations=speed_calculate(params);
create_plot('CMNC',x_dat(1:2),durations.CMNC,sel_ind)


disp('=====RANK COMPARISON ENDED =====')



%------------- M -------------
disp('===== VIEWS RANK COMPARISON STARTED=====')
clear all

size_all= {...
    {30,20,10},...
    {50,10},...
    {10,50}...
    };
noise_level = 0.01;
sparsity_level_all = [0.85];
k=1;
params.graph_tree(k) = graph_tree_root;
for i = 1:numel(size_all)
    params.graph_tree(k).Children(i).slices_num = 3;
    params.graph_tree(k).Children(i).noise_level = 0.01;
    params.graph_tree(k).Children(i).sparsity_level = 0.85;
    for j = 1:numel(size_all{i})
        params.graph_tree(k).Children(i).Children(j).type = 'clique';
        params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j}*2^1;
    end
end
params.graph_tree(k).labels = [];
[~,params.graph_tree(k)] = create_graph(params.graph_tree(k));
%     params.graph_tree(k)=cur_graph_tree;
% 
prms.embedding_method.GenClus.A_const = '+';
prms.embedding_method.GenClus.B_const = '+';
prms.embedding_method.GenClus.rho = 0;
prms.embedding_method.ComClus.beta = linspace(0.01,0.9145,3);
prms.embedding_method.ComClus.rho = linspace(0,0.16,3);
prms.embedding_method.ComClus.thres_inner = [1e-6];
prms.embedding_method.Symmetric_Richcom.structure = "true";
prms.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
prms.embedding_method.CMNC.structure = "true";
prms.embedding_method.CMNC.delta = 1;  

params.thres = [1e-6];
params.max_iters = [1000];
params.R = 3*2^5;
params.M = 3*2.^[0:5];
params.L_type_ind = [1 2];
params.sample = 1:5;


figure('Name','Rank (M)')
sel_ind = 4;
x_dat=params.M;

ax = axes;
hold on
ax.YScale='log';
ax.XScale='log';
ylabel('Duration (sec)')
xlabel('Views Rank (M)')
legend
grid

params.embedding_method=[];
params.embedding_method.GenClus=prms.embedding_method.GenClus;
durations=speed_calculate(params);
create_plot('GenClus',x_dat,durations.GenClus,sel_ind)

params.embedding_method=[];
params.embedding_method.ComClus=prms.embedding_method.ComClus;
durations=speed_calculate(params);
create_plot('ComClus',x_dat,durations.ComClus,sel_ind)

params.embedding_method=[];
params.embedding_method.Symmetric_Richcom=prms.embedding_method.Symmetric_Richcom;
durations=speed_calculate(params);
create_plot('Symmetric Richcom',x_dat,durations.Symmetric_Richcom,sel_ind)

params.embedding_method=[];
params.M=params.M(1:1);
params.embedding_method.CMNC=prms.embedding_method.CMNC;
durations=speed_calculate(params);
create_plot('CMNC',x_dat(1:1),durations.CMNC,sel_ind)


disp('===== VIEWS RANK COMPARISON ENDED =====')
end


function create_plot(method,x_dat,durations,sel_ind)
tmp = permute(durations,[sel_ind 1:(sel_ind-1) (sel_ind+1):numel(size(durations)) ]);
prctile_plot(x_dat,reshape(tmp,size(tmp,1),[])','DisplayName',method);
drawnow
end


function durations = speed_calculate(params)
embedding_method_print_type='nothing';
mtimesx_exists = exist('mtimesx','file');
alg_name_all=string(fieldnames(params.embedding_method));
max_iters=1000;
general_combs= [numel(alg_name_all) numel(params.graph_tree) numel(params.R) numel(params.M) numel(params.thres) numel(params.sample)];
if any(alg_name_all=="ComClus")
    comclus_combs= [general_combs(2:end)  numel(params.embedding_method.ComClus.beta) numel(params.embedding_method.ComClus.rho) numel(params.embedding_method.ComClus.thres_inner) ];
    comclus_combs_num=prod(comclus_combs);
end
comclus_ind=1;
comclus_total_time=0;
comclus_perc=0;
if any(alg_name_all=="Symmetric_Richcom")
    Symmetric_Richcom_combs= [general_combs(2:end)  numel(params.embedding_method.Symmetric_Richcom.rho) numel(params.embedding_method.Symmetric_Richcom.structure) ];
    Symmetric_Richcom_combs_num=prod(Symmetric_Richcom_combs);
end
Symmetric_Richcom_ind=1;
Symmetric_Richcom_total_time=0;
Symmetric_Richcom_perc=0;
if any(alg_name_all=="CMNC")
    CMNC_combs= [general_combs(2:end)  numel(params.embedding_method.CMNC.delta) numel(params.embedding_method.CMNC.structure) ];
    CMNC_combs_num=prod(CMNC_combs);
end
CMNC_ind=1;
CMNC_total_time=0;
CMNC_perc=0;
general_combs_num=prod(general_combs);
GenClus_combs_num=prod(general_combs(2:end));
GenClus_ind=1;
GenClus_total_time=0;
GenClus_perc=0;
total_time=0;
total_time_start=tic;
tmp='';
duration_embeddings_GenClus=[];
duration_embeddings_ComClus=[];
duration_embeddings_Symmetric_Richcom=[];
duration_embeddings_CMNC=[];
%                     fprintf('%d/%d %d/%d %d/%d %d/%d\n',graph_tree_ind,numel(params.graph_tree),R_ind,numel(params.R),M_ind,numel(params.M),thres_ind,numel(params.thres))
total_ind=0;
for ind = randperm(general_combs_num)
    [alg_name_ind, graph_tree_ind,R_ind,M_ind,thres_ind,sample]=ind2sub(general_combs,ind);
    graph_tree = params.graph_tree(graph_tree_ind);
    R = params.R(R_ind);
    M = params.M(M_ind);
    thres = params.thres(thres_ind);
    alg_name=alg_name_all(alg_name_ind);
    % GenClus
    GenClus_print_name='GenClus';
    CMNC_print_name='CMNC';
    ComClus_print_name='ComClus';
    Symmetric_Richcom_print_name='Richcom';
    switch alg_name
        case 'GenClus'
            GenClus_print_name = ['>' GenClus_print_name  '<'];
            alg=1;
            L_type_ind=2;
            cur_params=[]; 
            cur_params.embedding_method.GenClus.A_const = '+';
            cur_params.embedding_method.GenClus.B_const = '+';
            cur_params.embedding_method.GenClus.rho = 0;
            cur_embedding_method_name = string(fieldnames(cur_params.embedding_method));
            cur_alg_opts_names = fieldnames(cur_params.embedding_method.(cur_embedding_method_name));
            alg_opts = struct;
            for i=1:numel(cur_alg_opts_names)
                alg_opts.(cur_alg_opts_names{i})=cur_params.embedding_method.(cur_embedding_method_name).(cur_alg_opts_names{i});
            end
            [X,graph_tree] = create_graph(graph_tree);

            nodes_labels = {graph_tree.Children.labels};
            views_labels = graph_tree.labels;

            print_status()

            embeddings_real_time_start = tic;
            generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,embedding_method_print_type,mtimesx_exists);
            duration_embeddings_GenClus(sample,graph_tree_ind,R_ind,M_ind,thres_ind) = toc(embeddings_real_time_start);


            GenClus_total_time = GenClus_total_time + toc(embeddings_real_time_start);
            GenClus_perc = GenClus_ind/GenClus_combs_num*100;
            GenClus_ind = GenClus_ind +1;
        case 'ComClus'
            % ComClus
            ComClus_print_name = ['>' ComClus_print_name  '<'];
            alg=2;
            L_type_ind=1;
            cur_params=[]; 
            for beta_ind = randperm(numel(params.embedding_method.ComClus.beta))
                for rho_ind = randperm(numel(params.embedding_method.ComClus.rho))
                    for thres_inner_ind = randperm(numel(params.embedding_method.ComClus.thres_inner))
                        cur_params.embedding_method.ComClus.beta = params.embedding_method.ComClus.beta(beta_ind);
                        cur_params.embedding_method.ComClus.rho = params.embedding_method.ComClus.rho(rho_ind);
                        cur_params.embedding_method.ComClus.thres_inner = params.embedding_method.ComClus.thres_inner(thres_inner_ind);
                        cur_embedding_method_name = string(fieldnames(cur_params.embedding_method));
                        cur_alg_opts_names = fieldnames(cur_params.embedding_method.(cur_embedding_method_name));
                        alg_opts = struct;
                        for i=randperm(numel(cur_alg_opts_names))
                            alg_opts.(cur_alg_opts_names{i})=cur_params.embedding_method.(cur_embedding_method_name).(cur_alg_opts_names{i});
                        end
                        [X,graph_tree] = create_graph(graph_tree);
                        nodes_labels = {graph_tree.Children.labels};
                        views_labels = graph_tree.labels;

                        print_status()

                        embeddings_real_time_start = tic;
                        generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,embedding_method_print_type,mtimesx_exists);
                        duration_embeddings_ComClus(sample,graph_tree_ind,R_ind,M_ind,thres_ind,beta_ind,rho_ind,thres_inner_ind) = toc(embeddings_real_time_start);
                        comclus_total_time = comclus_total_time + toc(embeddings_real_time_start);
                        comclus_perc = comclus_ind/comclus_combs_num*100;
                        comclus_ind = comclus_ind +1;
                    end
                end
            end
        case 'Symmetric_Richcom'
            % Symmetric Richcom
            Symmetric_Richcom_print_name = ['>' Symmetric_Richcom_print_name  '<'];
            alg=3;
            L_type_ind=1;
            cur_params=[];
            for rho_ind = randperm(numel(params.embedding_method.Symmetric_Richcom.rho))
                for structure_ind = randperm(numel(params.embedding_method.Symmetric_Richcom.structure))
                    cur_params.embedding_method.Symmetric_Richcom.rho = params.embedding_method.Symmetric_Richcom.rho(rho_ind);
                    cur_params.embedding_method.Symmetric_Richcom.structure = params.embedding_method.Symmetric_Richcom.structure(structure_ind);
                    cur_embedding_method_name = string(fieldnames(cur_params.embedding_method));
                    cur_alg_opts_names = fieldnames(cur_params.embedding_method.(cur_embedding_method_name));
                    alg_opts = struct;
                    for i=randperm(numel(cur_alg_opts_names))
                        alg_opts.(cur_alg_opts_names{i})=cur_params.embedding_method.(cur_embedding_method_name).(cur_alg_opts_names{i});
                    end
                    [X,graph_tree] = create_graph(graph_tree);
                    nodes_labels = {graph_tree.Children.labels};
                    views_labels = graph_tree.labels;

                    print_status()

                    embeddings_real_time_start = tic;
                    generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,embedding_method_print_type,mtimesx_exists);
                    duration_embeddings_Symmetric_Richcom(sample,graph_tree_ind,R_ind,M_ind,thres_ind,rho_ind,structure_ind) = toc(embeddings_real_time_start);
                    Symmetric_Richcom_total_time = Symmetric_Richcom_total_time + toc(embeddings_real_time_start);
                    Symmetric_Richcom_perc = Symmetric_Richcom_ind/Symmetric_Richcom_combs_num*100;
                    Symmetric_Richcom_ind = Symmetric_Richcom_ind +1;
                end
            end
        case 'CMNC'
            % CMNC
            CMNC_print_name = ['>' CMNC_print_name  '<'];
            alg=4;
            L_type_ind=2;
            cur_params=[];
            for delta_ind = randperm(numel(params.embedding_method.CMNC.delta))
                for structure_ind = randperm(numel(params.embedding_method.CMNC.structure))
                    cur_params.embedding_method.CMNC.delta = params.embedding_method.CMNC.delta(delta_ind);
                    cur_params.embedding_method.CMNC.structure = params.embedding_method.CMNC.structure(structure_ind);
                    cur_embedding_method_name = string(fieldnames(cur_params.embedding_method));
                    cur_alg_opts_names = fieldnames(cur_params.embedding_method.(cur_embedding_method_name));
                    alg_opts = struct;
                    for i=randperm(numel(cur_alg_opts_names))
                        alg_opts.(cur_alg_opts_names{i})=cur_params.embedding_method.(cur_embedding_method_name).(cur_alg_opts_names{i});
                    end
                    [X,graph_tree] = create_graph(graph_tree);
                    nodes_labels = {graph_tree.Children.labels};
                    views_labels = graph_tree.labels;

                    print_status()

                    embeddings_real_time_start = tic;
                    generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,embedding_method_print_type,mtimesx_exists);
                    duration_embeddings_CMNC(sample,graph_tree_ind,R_ind,M_ind,thres_ind,delta_ind,structure_ind) = toc(embeddings_real_time_start);
                    CMNC_total_time = CMNC_total_time + toc(embeddings_real_time_start);
                    CMNC_perc = CMNC_ind/CMNC_combs_num*100;
                    CMNC_ind = CMNC_ind +1;
                end
            end
    end
    total_ind = total_ind + 1;
    print_status()
end

durations.GenClus = duration_embeddings_GenClus;
durations.ComClus = duration_embeddings_ComClus;
durations.Symmetric_Richcom = duration_embeddings_Symmetric_Richcom;
durations.CMNC = duration_embeddings_CMNC;
    function print_status()
        tmp2=tmp;
        tmp='';
        if any(alg_name_all=="GenClus")
            tmp = [tmp sprintf('%6.2f%%  %10s total time:  %s \t remaining time: \t %s\n',GenClus_perc,GenClus_print_name,duration(0,0,GenClus_total_time),duration(0,0,-GenClus_total_time*(1-1/GenClus_perc*100)))];
        end
        if any(alg_name_all=="ComClus")
            tmp =[tmp sprintf('%6.2f%%  %10s total time:  %s \t remaining time: \t %s\n',comclus_perc,ComClus_print_name,duration(0,0,comclus_total_time),duration(0,0,-comclus_total_time*(1-1/comclus_perc*100)))];
        end
        if any(alg_name_all=="Symmetric_Richcom")
            tmp =[tmp sprintf('%6.2f%%  %10s total time:  %s \t remaining time: \t %s\n',Symmetric_Richcom_perc,Symmetric_Richcom_print_name,duration(0,0,Symmetric_Richcom_total_time),duration(0,0,-Symmetric_Richcom_total_time*(1-1/Symmetric_Richcom_perc*100)))];
        end
        if any(alg_name_all=="CMNC")
            tmp =[tmp sprintf('%6.2f%%  %10s total time:  %s \t remaining time: \t %s\n',CMNC_perc,CMNC_print_name,duration(0,0,CMNC_total_time),duration(0,0,-CMNC_total_time*(1-1/CMNC_perc*100)))];
        end
        total_perc=total_ind/general_combs_num*100;
        tmp =[tmp sprintf('%6.2f%% \t    %10s:  %s\t remaining time: \t %s\n',...
            total_perc,'TOTAL TIME',duration(0,0,toc(total_time_start)),duration(0,0,-toc(total_time_start)*(1-1/total_perc*100)))];
    fprintf(repmat('\b',1,numel(tmp2)))
    fprintf("%s",tmp)
    end
end
