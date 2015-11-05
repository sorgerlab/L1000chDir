javaaddpath('mongo2.12.3.jar');

% wait for adding java path
pause(1);

import com.mongodb.MongoClient;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.util.JSON;
import java.util.Arrays;
import java.util.regex.Pattern;


mongoClient = MongoClient('10.91.53.225',27017);
db = mongoClient.getDB('LJP2014');

%
readColl = db.getCollection('experiments');
writeColl = db.getCollection('ljp56Chdirs');

% % test suit
% readColl = db.getCollection('ljp12Test');
% writeColl = db.getCollection('ljp12ChdirTest');

% utility function to convert java object to matlab struct.
j2m = @(x) loadjson(char(x.toString));


filter = BasicDBObject();

% % test suit
% filter.append('batch',Pattern.compile('LJP00[12]'));

filter.append('batch',Pattern.compile('LJP00[56]'));
batches = readColl.distinct('batch',filter);
batches = j2m(batches);


%% 

for i = 1%:numel(batches)
    batch = batches{i};
    fprintf('%s %d\n',batch,i);
    filter = BasicDBObject();
    filter.append('batch',batch);
    plates = readColl.distinct('det_plate',filter);
    plates = j2m(plates);
    
%     % replcate data on this plate has been removed from db
%     idx = strcmp('LJP005_MCF10A_3H_X2_B17',plates);
%     plates(idx) = [];
    platesRes = cell(numel(plates),1);
    
    % compute chdir for replicates on each plate
    for j = 1:numel(plates)
        plate = plates{j};
        query = BasicDBObject();
        query.append('det_plate',plate);
        cursor = readColl.find(query);
        arr = cell(cursor.count,1);
        for k = 1:cursor.count
            arr{k} = j2m(cursor.next());
        end
        % plateRes a struct: {arr:chdirsAndMeta,ctrlPMetricDistribution:vector}
        plateRes = getChdir(arr);
        platesRes{j} = plateRes;
    end
    
    % get unique experiments' sig_ids.
    jsonQuery = sprintf('[{$match:{batch:"%s",pert_type:{$ne:"ctl_vehicle"}}},{$group:{_id:{"batch":"$batch","pert_id":"$pert_id","pert_dose":"$pert_dose"},replicateCount:{$sum:1}}}]',batch);
    aggregateOutput = readColl.aggregate(JSON.parse(jsonQuery));
    sigIdStructs = j2m(aggregateOutput.results());
    
    % merge chdir replicates and computeSignficance
    chdirArr = mergeChdirReplicates(platesRes,sigIdStructs);
    
    % save chdir arr to db.
%     for j = 1:numel(chdirArr)
%         chdirStruct = chdirArr{j};
%         writeColl.save(BasicDBObject(JSON.parse(savejson('',chdirStruct))));
%     end
end

% notes: 
% add replicate count in the chdir field.
% use BasicDBObject(JSON.parse(jsonQuery)) for java query.
% sig_id for experiments is in the form "brew_prefix:pert_id:dose" (no dose unit).
% here batch == brew_prefix.
% distinct pert_type in LJP56: [ "trt_cp", "ctl_vehicle", "trt_poscon" ]


