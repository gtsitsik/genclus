function  clustering_quality_on_artificial_data()
% ---- Experiment setup ----

% Graphs with varying intra-community sparsity levels
size_all= {...
    {30,20,10},...
    {50,10},...
    {10,50}...
    };
noise_level = 0.01;
sparsity_level_all = [0.85:0.02:0.99];
for k = 1:numel(sparsity_level_all)
    params.graph_tree(k) = graph_tree_root;
    for i = 1:numel(size_all)
        params.graph_tree(k).Children(i).slices_num = 3;
        params.graph_tree(k).Children(i).noise_level = noise_level;
        params.graph_tree(k).Children(i).sparsity_level = sparsity_level_all(k);
        for j = 1:numel(size_all{i})
            params.graph_tree(k).Children(i).Children(j).type = 'clique';
            params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j}*2;
        end
    end
    params.graph_tree(k).labels = [];
    [~,params.graph_tree(k)] = create_graph(params.graph_tree(k));
%     params.graph_tree(k)=cur_graph_tree;
end

% Parameters that are common among all methods
params.thres = [1e-3 1e-6 1e-9];
params.max_iters = [1000];
params.R = 6:10;
params.M = 3;
params.L_type_ind = [1 2];
params.sample = 1:100;
params.workers = 20;
params.clustering_method.kmeans.replicates = [1];
params.clustering_method.kmeans.clusters_num = "3 2 2";
params.clustering_method.kmeans.row_normalization_type = ["none","unit"];  
params.clustering_method.large_inner_prod.thres = linspace(0.8,0.99,3);
params.clustering_method.maximum.nofield = "nofield";
params.column_normalization_type = ["none","sqrtB","B"];

% Sets up the parameters of GenClus, ComCLus and Symmetric Richcom.
params.embedding_method.GenClus.A_const = ['1','+','U'];
params.embedding_method.GenClus.B_const = ['1','+','U'];
params.embedding_method.GenClus.rho = 0;
params.embedding_method.ComClus.beta = linspace(0.01,0.9145,6);
params.embedding_method.ComClus.rho = linspace(0,0.16,6);
params.embedding_method.ComClus.thres_inner = [1e-6];
params.embedding_method.Symmetric_Richcom.structure = "true";
params.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
params.embedding_method.CMNC.structure = "true";
params.embedding_method.CMNC.delta = 1;

% Runs the experiment calculations
[~,~,~,~,~,filespath] = multi_run_calculate(params,true,true);


%---- Generation of figures ----

% Assigns a parameter to each visual dimension of the plots
xaxis = "graph_tree"; % horizontal axis
tile = ["clustering_measure", "clustered_entity"]; % horizontal and vertical plots
ln_clr = ["embedding_method"]; % line color
aggregate = "sample"; % shade corresponding to the 25-th and 75-th percentiles

% Plots for all enhanced methods 
fixed_params = [];
multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);
 
ln_clr = [""];

% Plots for original ComClus
fixed_params = ["clustering_method.choice = string('maximum')",
                "embedding_method.choice = string('ComClus')",
                "L_type_ind=1" ,
                "column_normalization_type=string('B')"];
multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);
 
% Plots for original Symmetric Richcom 
fixed_params = ["clustering_method.choice = string('large_inner_prod')",
                "embedding_method.choice = string('Symmetric_Richcom')",
                "L_type_ind=1",
                "column_normalization_type=string('none')"];
multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);
 
% Plots for original GenClus
fixed_params = ["embedding_method.choice = string('GenClus')",
                "L_type_ind=2",
                "embedding_method.GenClus.A_const = string('+')",
                "embedding_method.GenClus.B_const = string('+')"];
multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);

% Plots for original CMNC
fixed_params = ["clustering_method.choice = string('maximum')",
                "embedding_method.choice = string('CMNC')",
                "L_type_ind=2",
                "column_normalization_type=string('none')"];
multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);
