/* server.js
    Server initialization and server-side request handling code

    Requirements:   This code requires the node.js environment to function properly

    Jun Seok Lee <jlexternal@gmail.com> - May 2021
*/

/*
  Server-side requirements:
    An SQL database with the necessary tables. (see build_database.js)
    The modules imported below (and installed).

  Client-side requirements (if using Prolific):
    The URL to access the application MUST have a prolific_id query with
    a 24 character string specifying the Prolific ID.
*/

// import modules
var http        = require('http');
var fs          = require('fs');
var express     = require('express');
var path        = require('path');
var mysql       = require('mysql');
var bodyParser  = require('body-parser');

var buildDB     = require('./build_database');
const nsubj     = 100; // set to the number of CSV files to be entered into the db
const nblock    = 6;
const ntrial    = 80;

var queries_cache1      = [], // array to hold all queries to be printed to console (redundancy for safety)
    queries_cache2      = [],
    queries_cache_block = '',
    queries_cache_other = [];

var app = express();
app.set('port', process.env.PORT || 8000);                // listening port
app.use(express.static(path.join(__dirname, 'public')));  // adds directories to active path
app.use(bodyParser.urlencoded({ extended: true }));       // need to read JSON object body

console.log('New instance of server established.');

var db_config = {
  /* for use on local machine
  connectionLimit: 10,
  host: "localhost",
  user: "user",
  password: "password",
  database: 'database'
  //*/
};

// To use pooled SQL connections instead of individual ones
const pool = mysql.createPool(db_config);

// Promise-ified version of SQL query to get the task data
function sendQuery(conn,sql) {
  /*
    Note: This function is not in use. Kept for future reference.
  */
    return new Promise ( (resolve, reject) => {
      conn.query(sql, (error, results, fields) => {
        // SQL error handling
        if(error) {
          return reject(error);
        }
        resolve({results, fields});
      });
    });
}

function sendQueryPool(sql) {
  return new Promise ((resolve, reject) => {
    pool.getConnection(function(err, connection) {
      if (err) throw err; // not connected!

      connection.query(sql, (error, results, fields) => {
        connection.release(); // when done with the connection, release it.

        // SQL error handling
        if (error) {
          return reject(error);
        }
        resolve({results, fields});
      });
    });
  });
}

// write to database every 30 seconds
var sendCacheToggle     = 0;
const sendCacheInterval = 30000;
// query batch to database and print for redundancy/backup
function sendCachedQueries() {
  var queries_cache_str = 'INSERT INTO resp_table '+
                          '(unique_id, cond, i_block, i_round, i_trial, seen_feedback, is_correct, choice_position, choice_symbol, reaction_time) VALUES ';

  (new Promise((resolve,reject) => {
    // check the 1st cache
    if (sendCacheToggle == 0) {
      if (queries_cache1.length != 0) {
          //console.log('\n Sending queries from cache 1...'); //debug
          let query_str = queries_cache_str;
          sendCacheToggle = 1;
          queries_cache1.forEach(value => {
            query_str += value;
          });
          query_str = query_str.slice(0,-1) + ';';
          sendQueryPool(query_str)
          .then(resolve());
      } else {
          resolve();
      }
    }
    else {
      if (queries_cache2.length != 0) {
          //console.log('\n Sending queries from cache 2...'); //debug
          let query_str = queries_cache_str;
          sendCacheToggle = 0;
          queries_cache2.forEach(value => {
            query_str += value;
          });
          query_str = query_str.slice(0,-1) + ';';
          sendQueryPool(query_str)
          .then(resolve());
      } else {
        resolve();
      }
    }
  }))
  .then(()=> {
    // reset sent cache
    if (sendCacheToggle == 1) {
        queries_cache2 = [];
    } else {
        queries_cache1 = [];
    }
  })
  .catch((error) => console.log(error));
}
// send cached queries to db and then dump them
setInterval(function() {
  sendCachedQueries();
},sendCacheInterval);

// if launching the app for the first time, build the tables of the database
async function checkTablesExistence() {
  var check_tables_str = "SHOW TABLES LIKE 'main_table';";
  var qres = await sendQueryPool(check_tables_str);
  var ntables = 0;
  if (qres.results.length == 1) {
    ntables++;
    check_tables_str = "SHOW TABLES LIKE 'task_table';";
    qres = await sendQueryPool(check_tables_str);
    if (qres.results.length == 1) {
      ntables++;
      check_tables_str = "SHOW TABLES LIKE 'resp_table';";
      qres = await sendQueryPool(check_tables_str);
      if (qres.results.length == 1) {
        ntables++;
      }
    }
  }
  console.log('Tables checked');
  return ntables;
}

// check if database is built already, if not, build_database; if so, skip
checkTablesExistence()
.then(function(ntables) {
  if (ntables == 0) {
    console.log('No tables found. Building tables...');
    new Promise(resolve => {
      resolve(buildDB.create_tables());
    }).then(async(queryArray) => {
      //debug create tables upon launch
        var promises = [];
        for (var query of queryArray) {
          promises.push(sendQueryPool(query));
        }
        var isdone = await Promise.all(promises);
        return isdone;
    }).then(async() => {
        var main_table_queryArray = await buildDB.populate_main_table(nsubj);
        var task_table_queryArray = await buildDB.populate_task_table(nsubj,nblock,ntrial);
        var promises = [];
        for (var query1 of main_table_queryArray) {
          promises.push(sendQueryPool(query1));
        }
        for (var query2 of task_table_queryArray) {
          promises.push(sendQueryPool(query2));
        }
        var isdone = await Promise.all(promises);
    }).then(()=> {
        console.log('Connection ended (build_tables)');
        console.log('You should manually ensure that task_table is well filled.');
    });
  } else if (ntables == 1 | ntables == 2) {
      console.error('Serious error! There are some tables missing! Check database.');
      process.abort();
  } else {
      console.log('All tables found. Awaiting URL request...');
  }
});

// launch Express server
var server = http.createServer(app);
server.listen(app.get('port'), function () {
  console.log("Express server listening on port " + app.get('port'));
});

// create a route for the main page
app.get('/', (req, res) => {
  // get URL parameters (for Prolific ID)
  var pid = req.query.prolific_id; // The case of prolific_id is important. It should be in lowercase, or there may be problems with parsing it

  // make sure that the Prolific ID is exactly 24 characters
  if (pid.length != 24) {
    console.log('The queried Prolific ID is not exactly 24 characters in length! Aborting execution...');
  }

  // post IP address of participant (for debugging purposes)
  var ip = req.headers['x-forwarded-for'] || req.connection.remoteAddress;
  console.log('Client with prolific_id = '+pid+' has IP address ' + ip);

  // read the html file containing the experiment
  fs.readFile('public/run_expe_retest.html', function(err, data) {
    return res.end();
  });
  // send to client
  res.sendfile('public/run_expe_retest.html');
});

// SQL queries to update and get information at start of task
async function runGetTaskLoop(nb,nt,pid) {
  // check for repeat of prolific id
  let query_str_check = "SELECT EXISTS(SELECT * FROM main_table WHERE prolific_id = '"+pid+"');";

  let result = await sendQueryPool(query_str_check);

  Object.values(JSON.parse(JSON.stringify(result.results[0]))).map(val => {
    if (val != 0) { // if the EXISTS check is true
        console.log('Duplicate PID case!');
        console.log('Stopping operation for '+pid+'...');
        return false;
    } else {
        console.log('Moving on...');
    }
   });
   //debug: need to send an alert to the client saying that there is a repeat of the prolific id

  // get the first unique_id from server that is not already taken
  let query_str_uid = 'SELECT unique_id, isubj FROM main_table WHERE prolific_id IS NULL LIMIT 1;';
  let query_out_uid = await sendQueryPool(query_str_uid);
  let unique_id     = query_out_uid.results[0].unique_id;
  subj_num          = query_out_uid.results[0].isubj;

  // update main_table with Prolific ID and timestamp start of task
  let query_str_pid = 'UPDATE main_table ' +
                      "SET prolific_id = '" + pid + "', timestamp_task_start = NOW() " +
                      "WHERE unique_id = '" + unique_id + "';";
  try {
      await sendQueryPool(query_str_pid);
      console.log('Organizing task for UID '+ unique_id + '...');
  } catch(err) {
      console.error('Error in querying entry for '+pid+'... Cannot update main_table!');
      return false;
  }

  // get relevant task data for the specific participant
  var table_name = 'task_table';
  var traj_all = [];
  var idx_blocks = [];

  for (var ib=1; ib<nb+1; ib++) {
    let traj_block = [];
    let idx_block= [];
    for (var it=1; it<nt+1; it++) {
      // get feedback values for the correct option
      let query_str_traj = 'SELECT correct_choice_feedback FROM '+table_name+' WHERE ' +
                           "unique_id='"+unique_id+"' AND i_block="+ib+" AND i_trial="+it+";";
      let query_out_traj = await sendQueryPool(query_str_traj);
      traj_block.push(query_out_traj.results[0].correct_choice_feedback);

      // get rounds values for the 2 sessions (halves) of the experiment
      if (ib==1 | ib==4) {
        let query_str_idx = 'SELECT i_round FROM '+table_name+' WHERE ' +
                            "unique_id='"+unique_id+"' AND i_block="+ib+" AND i_trial="+it+";";

        let query_out_idx = await sendQueryPool(query_str_idx);
        idx_block.push(query_out_idx.results[0].i_round);
      }
    }
    traj_all.push(traj_block);
    if (ib==1 | ib==4) {
      idx_blocks.push(idx_block); // global idx_block filled here
    }
  }
  // create the task object to send to function
  var task_obj = {
    idx_blocks: idx_blocks,
    traj_all: traj_all,
    unique_id: unique_id,
    subj_num: subj_num
  };
  return task_obj;
}

// listen for POST request_task request from frontend
app.post('/request_task', (req, res) => {
  console.log('POST request_task requested from client!');
  const nblocks   = 6;
  const ntrials   = 80;

  if (req.body.pid) {
    let pid = req.body.pid;
    // variables to send to client
    let traj, idx_blocks, uid, subj_num;

    // loop through the blocks and trials to get task data
    runGetTaskLoop(nblocks,ntrials,pid)
      .then(function(result) {
        // error
        if (!result) {
          return false;
        } else {
          // fill variables to send
          traj       = result.traj_all;
          idx_blocks = result.idx_blocks;
          uid        = result.unique_id;
          subj_num   = result.subj_num;
          return true;
        }
      }) // global traj filled here
      .then(function(continueFlag) {
        if (continueFlag) {
          setTimeout(function() { // wait a little bit for the query to be processed
            var task_data = {
              traj: traj,
              idx_blocks: idx_blocks,
              unique_id: uid,
              subj_num: subj_num
            };
            res.send(task_data); // send id data to client
            res.end();
          },100);
        }
      })
      .catch(console.error); //debug : maybe return an erorr page to client
  }
});

// listen for server calls to write TASK RESPONSES to SQL database
app.post('/post_resp', (req, res) => {
  //console.log('POST resp requested from client!');
  if (req.body) {
    var table_name = 'resp_table'; // change to necessary table name
    var data = req.body,
     unique_id = data.unique_id, // CHAR
            cd = data.cond, // CHAR
            ib = data.i_block,
            ir = data.i_round,
            it = data.i_trial,
            fb = data.seen_feedback,
            cr = data.is_correct, //BOOL
            cp = data.choice_position, // CHAR
            cs = data.choice_symbol,
            rt = data.reaction_time;

    var char_fields = new Set(['unique_id','cond','choice_position']);
    var field_str = '';
    var value_str = '';
    for (var field of Object.keys(data)) {
      // concatenate list of table fields
      field_str += field+', ';
      // concatenate list of values to be input
      if (char_fields.has(field)) {
        // deal w/ CHAR
        value_str += "'" + data[field] + "', ";
      } else if (field == 'is_correct') {
        // deal w/ BOOL
        if (data[field] == 'true') {
          value_str += '1, ';
        } else {
          value_str += '0, ';
        }
      } else {
        // everything else is INT
        value_str += data[field] + ", ";
      }
    }
    value_str = '('+value_str.slice(0,-2)+'),';

    // queue queries to be batch sent
    if (sendCacheToggle == 0) {
        queries_cache1.push(value_str);
    } else {
        queries_cache2.push(value_str);
    }

    // print long insert query at end of each half-block
    if (it % 40 == 0) {
      queries_cache_block += value_str;
      (new Promise((resolve, reject) => {
        let insert_str = 'INSERT INTO resp_table '+
                         '(unique_id, cond, i_block, i_round, i_trial, seen_feedback, is_correct, choice_position, choice_symbol, reaction_time) VALUES ';
        console.log('Hemi-block response log dump:\n');
        console.log(insert_str.concat(queries_cache_block));
        console.log('\n');
      }))
      .then(queries_cache_block = '');
    } else {
        queries_cache_block += value_str;
    }
    res.send(true);
    res.end();
  } else {
      console.log('Parsed object body is empty!');
      res.send(false);
      res.end();
  }
});

// listen for server calls to write PARTICIPANT FEEDBACK to SQL database
app.post('/post_fb', (req, res) => {
  if (req.body) {
    var table_name = 'main_table';
    var data = req.body,
        unique_id = data.unique_id,
        feedback1 = data.entry1,
        feedback2 = mysql.escape(data.entry2);
    var feedback = feedback1.concat('Q0:"' +feedback2.slice(1,-1)+'"');
    var query_str = "UPDATE "+table_name+" SET feedback = '"+feedback+"' "+
                    "WHERE unique_id = '" + unique_id + "';";
    console.log('Prolific id: ' + unique_id + ' sent feedback!');
    queries_cache_other.push(query_str);
    sendQueryPool(query_str)
    .then(()=> {
      console.log('Connection released (post_fb)');
    });
    res.send(true);
  }
  res.end();
});

// update participant for bonuses if they achieve minimum acceptable score
app.post('/bonus', (req, res) => {
  if (req.body) {
    var unique_id = req.body.unique_id,
        query_str = "UPDATE main_table SET bonus_flag = 1 WHERE unique_id = '" + unique_id + "';";
    console.log('Prolific id: ' + unique_id + ' received a bonus!');
    console.log(query_str);
    queries_cache_other.push(query_str);
    sendQueryPool(query_str)
    .then(()=> {
      console.log('Connection released (bonus)');
    });
    res.send(true);
  }
  res.end();
});

// update database with flag and time when task ends
app.post('/end_task', (req, res) => {
  if (req.body) {
    var unique_id = req.body.unique_id;
    var query_str_end_task = 'UPDATE main_table ' +
                             "SET timestamp_task_end = NOW(), task_completion_flag = TRUE " +
                             "WHERE unique_id = '" + unique_id + "';";
    console.log('Prolific id: ' + unique_id + ' has ended the task!');
    console.log(query_str_end_task);
    queries_cache_other.push(query_str_end_task);
    sendQueryPool(query_str_end_task)
    .then(()=> {
      console.log('Connection released (end_task)');
    });
    res.send(true); // return 'true' back to client
  }
  res.end();
});

/*                                         */
/* server calls for the questionnaire page */
/*                                         */

// function to obtain id parameters from the SQL database (used during questionnaire)
async function getIDForQues(pid) {
  // get the unique_id and subj_num for the Prolific ID sent by client
  let query_str_uid   = "SELECT unique_id, isubj FROM main_table WHERE prolific_id = '" + pid +"';";
  let query_out_uid   = await sendQueryPool(query_str_uid);
  try {
      var idObj = {
        uid_ques: query_out_uid.results[0].unique_id,
        subj_num_ques: query_out_uid.results[0].isubj
      };
  } catch(err) {
      console.log(err);
  }
  return idObj;
}

// if participant has already finished the questionnaire, redirect to error
function checkAllowQues(pid) {
  return new Promise((resolve,reject) => {
    // check to see whether participant is eligible for accessing the questionnaire
    let query_str_check = "SELECT EXISTS(SELECT * FROM main_table WHERE prolific_id = '"+pid+"' AND timestamp_ques_end IS NULL);";
    try {
        sendQueryPool(query_str_check)
        .then((result) =>  {
          Object.values(JSON.parse(JSON.stringify(result.results[0]))).map(val => {
            if (val == 1) { // if the EXISTS check is true
                console.log('Allowing access to questionnaire for pid '+pid+' ...');
                resolve(true);
            } else {
                console.log('PID '+pid+' has already accessed the questionnaire!');
                resolve(false);
            }
          });
        });
    } catch (error) {
        return reject(error);
    }
  });
}

// create a route for the questionnaire page
app.get('/ques', (req, res) => {
  var pid = req.query.prolific_id;
  var ip = req.headers['x-forwarded-for'] || req.connection.remoteAddress;
  console.log('Client with prolific_id = '+pid+' has IP address ' + ip + 'is trying to access the questionnaire!');

  checkAllowQues(pid)
  .then((allowFlag) => {
    // make sure that the Prolific ID is exactly 24 characters
    if (pid.length != 24 || allowFlag == false) {
        if (pid.length != 24) {
          console.log('The queried Prolific ID is not exactly 24 characters in length!');
        }
        // read the html file containing the error
        fs.readFile('public/run_ques_error.html', function(err, data) {
          return res.end();
        });
        // send to client
        res.sendfile('public/run_ques_error.html');
    } else {
        // read the html file containing the questionnaire
        fs.readFile('public/run_ques.html', function(err, data) {
          return res.end();
        });
        // send to client
        res.sendfile('public/run_ques.html');
    }
  });
});

// gets id parameters at start of questionnaire
app.post('/request_id', (req,res) => {
  // if there's nothing we need to send this information to the client
  if (req.body.pid) {
    let pid = req.body.pid;
    getIDForQues(pid)
    .then(function(idObj) {
      // set questionnaire start timestamp
      let query_str_ques_start = 'UPDATE main_table SET timestamp_ques_start = NOW() '+
                                 'WHERE isubj = '+idObj.subj_num_ques + ';';
      sendQueryPool(query_str_ques_start);
      // wait a little bit for the query to be processed
      setTimeout(function() {
        res.send(idObj);
        res.end();
      },100);
    })
    .catch(console.error);
  }
});

// listen for server calls to write QUESTIONNAIRE RESPONSES to SQL database
app.post('/post_ques', (req, res) => {
  if (req.body) {
    var table_name = 'ques_table';
    var data = req.body,
      unique_id = data.unique_id,
      label     = data.label,
      responses = data.responses;
    var value_str = "'"+unique_id+"'"+', '+"'"+label+"'"+', '+"'"+responses+"'";
    var query_str = "INSERT INTO " +table_name+" (unique_id, label, responses) VALUES ("+value_str+");";
    console.log('Prolific id: ' + unique_id + ' posted their questionnaire!');
    console.log(query_str);
    queries_cache_other.push(query_str);
    sendQueryPool(query_str)
    .then(()=> {
      console.log('Connection released (post_ques)');
    });
    res.send(true);
  }
  res.end();
});

// update database with flag and time when questionnaire is completed
app.post('/end_questionnaire', (req, res) => {
  if (req.body) {
    var data = req.body;
    var unique_id = data.unique_id;
    var query_str_end_q = 'UPDATE main_table ' +
                          "SET timestamp_ques_end = NOW(), questionnaire_completion_flag = TRUE " +
                          "WHERE unique_id = '" + unique_id + "';";
    console.log('Prolific id: ' + unique_id + ' posted their questionnaire!');
    console.log(query_str_end_q);
    queries_cache_other.push(query_str_end_q);
    sendQueryPool(query_str_end_q)
    .then(()=> {
      console.log('Connection released (post_ques)');
    });
    res.send(true);
  }
  res.end();
});
