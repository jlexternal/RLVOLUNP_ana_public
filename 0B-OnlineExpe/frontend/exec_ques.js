/* exec_ques
    Runs the function executing the questionnaire.

    Requirements: qns - JSON object containing shuffled questionnaire battery (see qn_630.js)

    Jun Seok Lee <jlexternal@gmail.com> - May 2021
*/
function exec_ques(qns) {

  // present consent form
  var check_consent = function(elem) {
    if (document.getElementById('consent_checkbox').checked) {
      return true;
    }
    else {
      alert("If you wish to participate, you must check the box next to the statement 'I agree to participate in this study.'");
      return false;
    }
    return false;
  };
  // declare the consent block
  var consent_trial = {
    type:'external-html',
    url: "consent_form/consentpg_rlvolunp_ques.html",
    cont_btn: "startbutton",
    check_fn: check_consent
  };
  timeline.push(consent_trial);

  // push each questionnaire into the timeline
  qns.forEach(q => timeline.push(q)); //timeline is declared outside function

  // match the question in questionnaire to the scale for ordinal conversion
  function findScale(label,ques_num,iset) {
    var scale;
    switch (label) {
      case 'bis':
        scale = ["Do not agree at all",  "Agree slightly", "Agree a lot",  "Agree completely"];
        break;
      case 'ocir':
        scale = ["Not at all", "A little", "Moderately", "A lot", "Extremely"];
        break;
      case 'schizo':
        scale = ["No", "Yes"];
        break;
      case 'depress':
        scale = ["A little of the time",  "Some of the time", "Good part of the time",  "Most of the time"];
        break;
      case 'social':
        if (ques_num % 2 == 0) {
            scale = [ "None", "Mild", "Moderate", "Severe"];
        } else {
            scale = [ "Never (0%)", "Occasionally (1\-33%)", "Often (33\-67%)", "Usually (67\-100%)"];
        }
        break;
      case 'iq':
        var iq_scales = {
          s1:  ["2", "3", "4", "5", "6","7"],
          s2:  ["Richard is taller than Matt", "Richard is shorter than Matt",  "Richard is as tall as Matt" ,"It's impossible to tell"],
          s3:  ["25" ,"39", "44", "47" ,"53", "57"],
          s4:  ["Friday" ,"Monday", "Wednesday",  "Saturday" ,"Tuesday", "Sunday"],
          s5:  ["S" ,"T", "U",  "V" ,"W", "X"],
          s6:  ["E" ,"F", "G", "H" ,"I", "J"],
          s7:  ["T" ,"U", "V",  "X" ,"Y", "Z"],
          s8:  ["J" ,"H", "I", "N" ,"M", "L"],
          s10: ["A" ,"B", "C", "D" , "E", "F","G","H"]
        };
        if (iset == 1) {
            scale = iq_scales.s10;
        } else {
            scale = iq_scales['s'+(ques_num+1)];
        }
        break;
      case 'alcohol':
        var alc_scales = {
          s1: [ "Never" , "Monthly or less", "Two to four times a month","Two to three times a week","Four or more times a week"],
          s2: [ "1 or 2", "3 or 4", "5 or 6", "7 to 9", "10 or more"],
          s3: [ "Never", "Less than monthly", "Monthly", "Weekly", "Daily or almost daily"],
          s4: [ "No", "Yes, but not in the last year", "Yes, during the last year"]
        };
        if (ques_num == 0) {
            scale = alc_scales.s1;
        } else if (ques_num == 1) {
            scale = alc_scales.s2;
        } else if (ques_num > 1 & ques_num < 8) {
            scale = alc_scales.s3;
        } else {
            scale = alc_scales.s4;
        }
        break;
      case 'apathy':
        scale = [ "Not at all characteristic", "Slightly characteristic", "Somewhat characteristic", "Very characteristic"];
        break;
      case 'eat':
        scale = [ "Always", "Usually", "Often", "Sometimes", "Rarely", "Never"];
        break;
      case 'anxiety':
        scale = ["Almost never", "Sometimes", "Often", "Almost always"];
        break;
    }
    return scale;
  } // findScale

  // converts a string response to an ordinal (Likert scale, starting from 1)
  function convertToOrdinal(resp,scale) {
    var ord_resp = scale.indexOf(resp);
    return ord_resp + 1;
  }

  function post_ques(data) {
    var label   = data.label;
    // ignore welcome trial
    if (typeof label == 'undefined') {
      return;
    }
    // create object to send to server
    var sendObj = {
      unique_id: unique_id, // unique_id assigned outside of function
      label: label,
      responses: ''
    };
    // process response data and convert
    for (i=0; i<Object.keys(data.response).length; i++) {
      var iset = 0;
      if (label == 'iq') {
        const iqRespSet2 = new Set(["A" ,"B", "C", "D" , "E", "F","G","H"]);
        if (iqRespSet2.has(Object.values(data.response)[0])) {
          iset = 1;
        }
      }
      let scale = findScale(label,i,iset);
      let ord_resp = convertToOrdinal(data.response['Q'+i],scale); // convert to Likert-scaled response
      sendObj.responses = sendObj.responses.concat('"Q'+i+':"'+ord_resp+'",'); // write out the response string
    }
    // send response data to server
    $.ajax({
      url: '/post_ques',
      type: 'POST',
      data : sendObj,
  //  success: function(response){console.log('SUCCESS response POST request!');},  //debug use
      error: function(error){
        console.log(error);
      }
    })
  //  .done((successFlag) => { console.log('response POST successFlag: '+ successFlag); }) //debug use
      .fail((xhr,txt,err) => { console.log('FAIL post_ques POST'); console.log(err); });
  }

  // end of questionnaire prompt
  var end_questionnaire = {
    type: 'html-keyboard-response',
    stimulus: 'End of questionnaire. <br><br>'+
              'Thank you for your participation.<br><br>'+
              'You will be redirected to Prolific in a few moments... <span style = "font-weight:bold">Do not close the browser!</span>',
    choices: jsPsych.NO_KEYS,
    trial_duration: 3000
  };
  timeline.push(end_questionnaire);

  var labels = new Set(['bis', 'ocir', 'schizo', 'depress', 'social', 'iq', 'alcohol', 'apathy', 'eat', 'anxiety']);

  // jsPsych general handling parameters
  jsPsych.init({
    timeline: timeline,
    minimum_valid_rt: 100,
    on_trial_finish: function(data) {
      // SQL database query here:
      if (labels.has(data.label)) {
        post_ques(data);
      }
    },
    on_finish: function() {
      // send AJAX query to update main_table with questionnaire_completion_flag to TRUE
      $.ajax({
        url: '/end_questionnaire',
        type: 'POST',
        data: {unique_id: unique_id},
    //  success: function(response){console.log('SUCCESS response POST request!');},  //debug use
        error: function(error){
          console.log(error);
        }
      })
        .done((successFlag) => {
            window.location = "https://app.prolific.co/submissions/complete?cc=8F7093B9"; //debug : set to prolific completion page for QUESTIONNAIRE
        })
        .fail((xhr,txt,err) => { console.log('FAIL post_ques POST'); console.log(err); });
    }
  });

}
