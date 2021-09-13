# 2-OnlineExpe
Code for executing the online task.

**Note: This code will not run right out the gate. It requires a MySQL database and the proper communication channels 
to be opened between it and the server.** 

I personally used a Node.js server and an Express application framework as indicated in server.js. 
See https://nodejs.org/ for more information.

## Before running the experiment
Running the online experiment requires 3 parts, the backend, the frontend, and the database.

### The database
Given that you have a working MySQL database communicating with the server, and if all the task data is 
in the folder /public/csvs/ and the number of tasks are matched inside server.js, the tables for your database
will automatically be filled upon first execution of the server.

It is not necessary to use MySQL, but then you will have to write and setup your own architecture. 
See the file build_database.js for the details concerning the contents of the database.

### The server
As mentioned above, the server I used was Node.js, but you may use whatever as long as it can communicate with 
the database and the frontend. 
All details concerning the server can be found in the file server.js.

### The frontend
The frontend is the actual experiment. With a working backend server and database, the frontend will execute without modification.
It uses AJAX queries to communicate with the server, and these are not particular to Node.js.
