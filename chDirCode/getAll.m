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


mongoClient = MongoClient('localhost',27017);
db = mongoClient.getDB('LJP2014');

%
readColl = db.getCollection('ljp56Chdirs');

% utility function to convert java object to matlab struct.
j2m = @(x) loadjson(char(x.toString));
%%

cursor = readColl.find();
ljp56ChdirAll = cell(cursor.count(),1);
for i = 1:cursor.count()
    if mod(i,50) == 0
        disp(i);
    end
   ljp56ChdirAll{i} = j2m(cursor.next());
end
