function real_world_case_study()

% Parse input data
routes = readtable('routes.dat', 'ReadVariableNames', false);

routes.Properties.VariableNames = {'Airline','Airline ID','Source airport','Source airport ID','Destination airport','Destination airport ID','Codeshare','Stops','Equipment'};
routes.("Source airport ID") = string(routes.("Source airport ID"));
routes.("Destination airport ID") = string(routes.("Destination airport ID"));
routes.("Airline ID") = string(routes.("Airline ID"));
airports =  readtable('airports.dat', 'ReadVariableNames', false,'Format','auto');
airports.Properties.VariableNames = {'Airport ID','Name','City','Country','IATA','ICAO','Latitude','Longitude','Altitude','Timezone','DST','Tz database time zone','Type','Source'};
airports.("Airport ID") = string(airports.("Airport ID"));
airlines =  readtable('airlines.dat', 'ReadVariableNames', false);
airlines.Properties.VariableNames = {'Airline ID','Name','Alias','IATA','ICAO','Callsign','Country','Active'};
airlines.("Airline ID")= string(airlines.("Airline ID"));
countries =  readtable('countries.dat', 'ReadVariableNames', false);
countries.Properties.VariableNames = {'name','iso_code','dafif_code'};

countries.dafif_code = [];
countries = unique(countries);

continents = readtable('continents2.csv') ;
airlines_ext = airlines(logical(prod(string(airlines.Country)~=string(setdiff(airlines.Country,countries.name))',2)),:);
airlines_ext2 = join(airlines_ext,countries(:,{'name','iso_code'}),'LeftKeys','Country','RightKeys','name');
airlines_ext3 = airlines_ext2(logical(prod(string(airlines_ext2.iso_code)~=string(setdiff(airlines_ext2.iso_code,continents.alpha_2))',2)),:);
airlines_ext4 = join(airlines_ext3,continents(:,{'alpha_2','region'}),'LeftKeys','iso_code','RightKeys','alpha_2');

airports_ext = airports(logical(prod(string(airports.Country)~=string(setdiff(airports.Country,countries.name))',2)),:);
airports_ext2 = join(airports_ext,countries(:,{'name','iso_code'}),'LeftKeys','Country','RightKeys','name');
airports_ext3 = airports_ext2(logical(prod(string(airports_ext2.iso_code)~=string(setdiff(airports_ext2.iso_code,continents.alpha_2))',2)),:);
airports_ext4 = join(airports_ext3,continents(:,{'alpha_2','region'}),'LeftKeys','iso_code','RightKeys','alpha_2');

routes(routes.("Source airport ID")=="\N" | routes.("Destination airport ID")=="\N",:)=[];
tmp = union(setdiff(routes.("Source airport ID"),airports_ext4.("Airport ID")),setdiff(routes.("Destination airport ID"),airports_ext4.("Airport ID")));
routes = routes(all(routes.("Source airport ID")~=tmp',2) &all(routes.("Destination airport ID")~=tmp',2),:);
tmp = setdiff(routes.("Airline ID"),airlines_ext4.("Airline ID"));
routes = routes(all(routes.("Airline ID")~=tmp',2),:);

routes = join(routes,airports_ext4(:,{'Airport ID','region'}),'LeftKeys','Source airport ID','RightKeys','Airport ID');
routes.Properties.VariableNames{10} = 'Source Airport region';
routes = join(routes,airports_ext4(:,{'Airport ID','region'}),'LeftKeys','Destination airport ID','RightKeys','Airport ID');
routes.Properties.VariableNames{11} = 'Destination Airport region';
routes = join(routes,airlines_ext4(:,{'Airline ID','region'}),'LeftKeys','Airline ID','RightKeys','Airline ID');
routes.Properties.VariableNames{12} = 'Airline region';

european_airlines_routes = routes(logical(sum(routes.("Airline ID")==string(airlines_ext4(any(airlines_ext4.region==["Europe","Americas","Asia"],2),:).("Airline ID"))',2)),:);
european_airlines_routes = european_airlines_routes(logical(sum(european_airlines_routes.("Source airport ID")==string(airports_ext4(any(airports_ext4.region==["Europe","Americas","Asia"],2),:).("Airport ID"))',2)),:);
european_airlines_routes = european_airlines_routes(logical(sum(european_airlines_routes.("Destination airport ID")==string(airports_ext4(any(airports_ext4.region==["Europe","Americas","Asia"],2),:).("Airport ID"))',2)),:);

european_airlines_airport_IDs = table(string(intersect(european_airlines_routes.("Source airport ID"),european_airlines_routes.("Destination airport ID"))));
european_airlines_airport_IDs.Properties.VariableNames{1} = 'Airport ID';
european_airlines_airport_IDs = join(european_airlines_airport_IDs,airports_ext4(:,{'Airport ID','region'}),'LeftKeys','Airport ID','RightKeys','Airport ID');
european_airlines_airport_IDs = sortrows(european_airlines_airport_IDs,"region").("Airport ID");

european_airlines_IDs = table(unique(european_airlines_routes.("Airline ID")));
european_airlines_IDs.Properties.VariableNames{1} = 'Airline ID';
european_airlines_IDs = join(european_airlines_IDs,airlines_ext4(:,{'Airline ID','region'}),'LeftKeys','Airline ID','RightKeys','Airline ID');
european_airlines_IDs = sortrows(european_airlines_IDs,"region").("Airline ID");


% Create adjacency tensor
X = sptensor([],[],[numel(european_airlines_airport_IDs),numel(european_airlines_airport_IDs),numel(european_airlines_IDs)]);
for z = 1:size(european_airlines_routes,1)
    if mod(z,1000)==0
        disp(z+"/"+size(european_airlines_routes,1))
    end
    i = find(european_airlines_routes(z,:).("Source airport ID") == european_airlines_airport_IDs);
    j = find(european_airlines_routes(z,:).("Destination airport ID") == european_airlines_airport_IDs);
    k = find(european_airlines_routes(z,:).("Airline ID")==european_airlines_IDs);

    if ~isempty(i) &&  ~isempty(j) && ~isempty(k)
        X(i,j,k)= X(i,j,k)+1;
    end
end
valid_airport_inds=0;
valid_airline_inds=0;
european_airlines_airport_IDs2 =european_airlines_airport_IDs;
european_airlines_IDs2 =european_airlines_IDs;
X2=X;


% Data preprocessing as described in the appendix of the paper
while ~all(valid_airport_inds) || ~all(valid_airline_inds )
    if isa(X2,'double')
        valid_airline_inds = permute(sum(X2>0,[1,2]),[3,1,2])>100;
        valid_airport_inds = sum(X2>0,[2,3])>30;
    else
        valid_airline_inds=logical([]);
        valid_airport_inds=logical([]);
        for i =1:size(X2,3)
            valid_airline_inds(i)= nnz(X2(:,:,i))>100;
        end
        for i =1:size(X2,1)
            valid_airport_inds(i)= nnz(X2(i,:,:))>30;
        end
    end
    european_airlines_airport_IDs2 = european_airlines_airport_IDs2(valid_airport_inds);
    european_airlines_IDs2 = european_airlines_IDs2(valid_airline_inds);
    X2 = X2(find(valid_airport_inds),find(valid_airport_inds),find(valid_airline_inds));
    size(X2)
end
X2 = double(X2);

% Extract the labels of airlines and airports 
airlines_final = join(table(european_airlines_IDs2),airlines_ext4(:,{'Airline ID','region'}),'LeftKeys','european_airlines_IDs2','RightKeys','Airline ID');
airports_final = join(table(european_airlines_airport_IDs2),airports_ext4(:,{'Airport ID','region'}),'LeftKeys','european_airlines_airport_IDs2','RightKeys','Airport ID');


% Create a graph object suitable for Multi-Graph Explorer
par_all = graph_tree_root;
par_all.Data = X2;
tmp2 = string(airlines_final.region);
par_all.labels = tmp2; 
tmp2 = string(airports_final.region);
for i =1:numel(unique(tmp2))
    par_all.Children(i).labels = tmp2; 
end

% Run Multi-Graph Explorer
single_run_GUI(par_all)
