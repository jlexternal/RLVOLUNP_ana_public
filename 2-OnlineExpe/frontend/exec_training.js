/* exec_training
    Runs the function executing the training for each condition.

    Jun Seok Lee <jlexternal@gmail.com> - May 2021
*/
function exec_training(cond_type) {

  // the argument 'cond_type' determines which training session to run

  /*
    Introduce how the experiment works and what the goal of the task is.
    i.e. Not to choose the highest number everytime, but rather to choose
    the option that on average gives the higher value.

    Run a training block where the good option starts with a lower value but is the correct shape.

    Run a training block

    Familiarise the subject with the dynamics of the VOL condition as they will not be able to
    know where the switches happen.
  */

  /* load and preload images to be used as training stimuli */
  var shapes = ['img/shape_train01.png','img/shape_train02.png','img/shape_blank.png','img/shape_spacer.png'];
  var preload = {
    type: 'preload',
    images: shapes
  };
  timeline.push(preload);
  var stimulis = shapes.slice(0,2);

  // hard code reward values for the different types of training blocks
  var feedback_all = [
    /* REF */
    [[62, 59, 56, 48, 59, 73, 45, 54, 70, 67]],
    /* VOL */
    /*
      Note 1: The switch point is hard coded into the function debriefTableCreate.
              It will not be difficult to make this more dynamic with an argument.
      Note 2: All values in the array below are for the CORRECT choice

    */
    [[66, 41, 69, 50, 56, 55, 75, 78, /*switch point*/ 48, 67, 55, 50, 51, 77, 27, 54, 56, 63, 64, 74]],
    /* UNP */
     // values which fluctuate more wildly and hit closer to the extremes
    [[61, 62, 33, 57, 32, 46, 53, 67, 36, 68, 56, 54, 48, 41, 68, 57, 59, 72, 59, 43]]
  ];

  /* Load feedback for the respective training session */
  var fb_session;
  var cond_num;
  switch(cond_type) {
    case 'REF':
      fb_session = feedback_all[0];
      cond_num = 1;
      break;
    case 'VOL':
      fb_session = feedback_all[1];
      cond_num = 2;
      break;
    case 'UNP':
      fb_session = feedback_all[2];
      cond_num = 3;
      break;
  }

  var choice_opts = ['f','j'];          // set of available keys to choices
  var n_blocks    = fb_session.length;  // number of blocks in specified training session

  // preblock label
  function pre_training_block_label_fn(cond_num,iround,nround) {
    var pre_block_label = {
      // display round number if REF or UNP
      type: 'html-keyboard-response-faded',
      stimulus: function() {
        var stim_text;
        if (cond_num != 2) { // REF and UNP
          stim_text = '<p style = "font-size: 28px; font-weight: bold">Training</p><br><br>' +
          'Round ' + iround + '/' + nround;
        } else { // VOL
          stim_text = '<p style = "font-size: 28px; font-weight: bold">Training</p><br><br>' +
          'Single long round';
        }
        return stim_text;
      },
      minimum_duration: duration_fn(debugFlag,1000,100),
      choices: jsPsych.NO_KEYS,
      fadein_duration: duration_fn(debugFlag,2000,100),
      fadeout_duration: duration_fn(debugFlag,2000,100),
      trial_duration: duration_fn(debugFlag,1000,100),
    };
    return pre_block_label;
  }

  // block loop
  for (iblk=0; iblk<n_blocks; iblk++) {
    let fb_vals   = fb_session[iblk]; // array of feedback for the given block
    let n_trials  = fb_vals.length;   // number of trials for given block
    let curr_iblk = iblk;

    timeline.push(pre_training_block_label_fn(cond_num,iblk+1,n_blocks));

    /* Generate random placement of the correct choice
        The correct choice however, is always the shape_train01 */
    var cor_locs  = Array.from({length: n_trials}, () => Math.round(Math.random()));

    let score = 0;
    let score_arr = [];

    // trial loop
    for (itrl=0; itrl<n_trials; itrl++) {
      // localize indexed variables
      let stim_loc    = cor_locs[itrl];
      let cor_choice;
      if (cond_type == 'VOL' & itrl > 7) {
        cor_choice  = choice_opts[-stim_loc+1];
      } else {
        cor_choice  = choice_opts[stim_loc];
      }

      let stims       = [stimulis[stim_loc], stimulis[-stim_loc+1]];
      // main stimuli for each trial
      var trialstim = {
        type: "html-keyboard-response",
        stimulus: function(){
          return '<img src="' + stims[0]  + '"/><img src="img/shape_spacer.png" /><img src="' + stims[1]  + '" />';
        },
        choices: ['f','j'], // response set
        data: {
          task: 'response',
          correct_response: jsPsych.timelineVariable('correct_response'),
          condition: cond_type + '_training',
          trial_expe: itrl
        },
        on_finish: function (data){
          // correct/incorrect flag for trial
          if(jsPsych.pluginAPI.compareKeys(data.response, cor_choice)){
            data.correct = true;
            score++;
            score_arr.push(1);
          } else {
            data.correct = false;
            score_arr.push(0);
          }
        }
      };

      // update stimuli based on chosen option
      var chosen_trialstim = {
        type: "html-keyboard-response",
        stimulus: function(){
          if (jsPsych.data.get().last(1).values()[0].response == 'f') {
            return '<img src="' + stims[0]  + '"/><img src="img/shape_spacer.png" /><img src="img/shape_blank.png" />';
          }
          else {
            return '<img src="img/shape_blank.png"/><img src="img/shape_spacer.png" /><img src="' + stims[1]  + '" />';
          }
        },
        trial_duration: duration_fn(debugFlag,500,100),
        choices: jsPsych.NO_KEYS
      };

      // feedback stimulus
      let val = fb_vals[itrl];
      var feedback = {
        type: 'html-keyboard-response',
        stimulus: function(){
          /* *************************************************************************
          Warning: the argument of 'last()' must be the number of trials
            'feedback' is away from 'trialstim' in the timeline push below.
            e.g. 'timeline.push(trialstim, feedback)' ->  'last(1)'                     */

          let last_trial_correct = jsPsych.data.get().last(2).values()[0].correct;
          /* *************************************************************************  */
          let fb_str;
          if (last_trial_correct) {
            fb_str = val;
          } else {
            fb_str = 100-val;
          }
          return '<div style="font-size:60px;">'+fb_str+'</div>';
        },
        choices: jsPsych.NO_KEYS,
        trial_duration: duration_fn(debugFlag,1000,100),
      };

      // push events to timeline
      timeline.push(trialstim,chosen_trialstim,feedback);

      // end of last trial
      if (itrl == n_trials-1){

        var trialstim_end_feedback_total = {
          type: 'html-keyboard-response-faded',
          stimulus: function() {
            var fb_seen = [];
            for (var i=0; i<fb_session[curr_iblk].length; i++) {
              if (score_arr[i] == 1) {
                fb_seen.push(fb_session[curr_iblk][i]);
              } else {
                fb_seen.push(100-fb_session[curr_iblk][i]);
              }
            }
            // preface to the feedback table
            function boldANodeFn () {
              var boldASpanNode = document.createElement("span");
                  boldASpanNode.style.fontWeight = "bold";
                  boldASpanNode.appendChild(document.createTextNode("A"));
              return boldASpanNode;
            }
            var boldBSpanNode = document.createElement("span");
                boldBSpanNode.style.fontWeight = "bold";
                boldBSpanNode.appendChild(document.createTextNode("B"));

            let aboveTableSpanNode = document.createElement("span");
                aboveTableSpanNode.appendChild(document.createTextNode("Here's your performance during this training session:"));
                aboveTableSpanNode.appendChild(document.createElement("br"));
                aboveTableSpanNode.appendChild(document.createElement("br"));
                  greenTextSpanNode = document.createElement("span");
                  greenTextSpanNode.style.color = "green";
                  greenTextSpanNode.style.fontWeight  = "bold";
                  greenTextSpanNode.appendChild(document.createTextNode('O'));
                aboveTableSpanNode.appendChild(document.createTextNode("The green "));
                aboveTableSpanNode.appendChild(greenTextSpanNode);
                aboveTableSpanNode.appendChild(document.createTextNode(" indicates where you chose the better shape and "));
                aboveTableSpanNode.appendChild(document.createElement("br"));
                  redTextSpanNode = document.createElement("span");
                  redTextSpanNode.style.color = "red";
                  redTextSpanNode.style.fontWeight  = "bold";
                  redTextSpanNode.appendChild(document.createTextNode('X'));
                aboveTableSpanNode.appendChild(document.createTextNode("The red "));
                aboveTableSpanNode.appendChild(redTextSpanNode);
                aboveTableSpanNode.appendChild(document.createTextNode(" indicates where you chose the worse shape."));
                aboveTableSpanNode.appendChild(document.createElement("br"));
                aboveTableSpanNode.appendChild(document.createElement("br"));
                if (cond_type=='REF') {
                  aboveTableSpanNode.appendChild(document.createTextNode("As you can see, even if "));
                  aboveTableSpanNode.appendChild(boldBSpanNode); // B
                  aboveTableSpanNode.appendChild(document.createTextNode(" sometimes gave points higher than 50, or "));
                  aboveTableSpanNode.appendChild(boldANodeFn()); // A
                  aboveTableSpanNode.appendChild(document.createTextNode(" gave points lower than 50,"));
                  aboveTableSpanNode.appendChild(document.createElement("br"));
                  aboveTableSpanNode.appendChild(document.createTextNode("the good choice was shape "));
                  aboveTableSpanNode.appendChild(boldANodeFn()); // A
                  aboveTableSpanNode.appendChild(document.createTextNode("."));
                } else if (cond_type=='VOL') {

                } else { // UNP

                }
                aboveTableSpanNode.appendChild(document.createElement("br"));

            // create debriefing feedback table
            let condtype  = cond_type;
            let scorearr  = score_arr;
            let ntrials   = n_trials;
            let fbseen    = fb_seen;
            let tbl = debriefTableCreate(condtype,scorearr,ntrials,fbseen); // call the detailed feedback history function from above

            // post-feedback table commments
            let belowTableSpanNode = document.createElement("span");
            if (cond_type=='VOL') {
              belowTableSpanNode.appendChild(document.createTextNode("As you see above, the good choice switched to the other shape once during the round."));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createTextNode("During the real game, this will happen multiple times."));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createTextNode("Note, however, it will not flip over too fast after any given switch."));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createElement("br"));

            } else if (cond_type=='REF') {
              /* nothing */

            } else { // UNP
              belowTableSpanNode.appendChild(document.createTextNode("As you see, it may be more difficult to determine the good shape."));
              belowTableSpanNode.appendChild(document.createElement("br"));
              belowTableSpanNode.appendChild(document.createTextNode("Your gain is still based on the number of times you choose good shape."));
              belowTableSpanNode.appendChild(document.createElement("br"));
            }
            belowTableSpanNode.appendChild(document.createTextNode("Your theoretical true/paid score is: " + score + '/' + n_trials ));
            belowTableSpanNode.appendChild(document.createElement("br"));
            belowTableSpanNode.appendChild(document.createElement("br"));

            // continue instructions
            let continueSpanNode = document.createElement("span");
                continueSpanNode.style.fontWeight = 'bold';
                continueSpanNode.appendChild(document.createTextNode('Press spacebar to continue'));

            // stack the visuals to be displayed on screen
            let divDisplay = document.createElement("div"); //create new <div>
                divDisplay.appendChild(aboveTableSpanNode);
                divDisplay.appendChild(document.createElement("br"));
                divDisplay.appendChild(tbl);
                divDisplay.appendChild(document.createElement("br"));
                divDisplay.appendChild(document.createElement("br"));
                divDisplay.appendChild(belowTableSpanNode);
                divDisplay.appendChild(continueSpanNode);

            return divDisplay.innerHTML;
          },
          choices: [' '],
          fadein_duration: 500
        };
        timeline.push(trialstim_end_feedback_total);
      }

    } // end trial loop
  } // end block loop
}
