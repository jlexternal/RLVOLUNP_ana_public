function spit_long_instructions (instr_type) {

var instructions_var;
let instr_gen1, instr_gen2, instr_gen3, instr_gen4, instr_gen5, instr_gen6;

switch (instr_type) {
  case 'introduction':
    instr_gen1 = '<p style = "text-align: center; font-size: 28px">' +
      "To participate in the experiment you will mainly use two keys, <br><br>" +
      '<span style="font-weight:bold">F</span> and <span style="font-weight:bold">J</span> on your keyboard.</p><br>'+
      '<p style = "font-size: 24px; font-weight: bold">Press J to continue.</p>';
    instr_gen2 = '<p style = "text-align: center; font-size: 28px">' +
      "If you made it to this page, you've successfully clicked on "+'<span style = "font-weight: bold">J</span>!<br><br>' +
      '<p style = "font-size: 24px; font-weight: bold">Press F to go back.<br><br>Press J to continue.</p>';
    instr_gen3 = '<p style = "text-align: center; font-size: 28px; font-weight: bold">' +
      "Please read the following instructions very carefully.</p><br><br>" +
      '<span style="font-style:italic">Note: At the end of the experiment, you will be redirected to Prolific.<br>'+
      "If you leave the experiment at any point before this, payment cannot be guaranteed.</span><br><br>"+
      '<p style = "font-size: 24px; font-weight: bold">Press F to go back.<br><br>Press J to continue.</p>';
    instructions_var = {
      type: 'instructions-faded',
      pages: [instr_gen1, instr_gen2, instr_gen3],
      key_backward: 'f',
      key_forward: 'j',
      fade_duration: duration_fn(debugFlag,500,100),
    };
    break;

  case 'keypress':
    instr_gen1 = { stimuli: [
      'In our experiment, you will play a game where you choose from one of two choices represented by two symbols.<br><br>',
      "For now, we will represent the two symbols with letters: <br>"+
      '<span style="font-weight:bold">🄰</span> and '+'<span style="font-weight:bold">🄱.</span><br><br>',
      // The letters in boxes is Unicode text, if it shows up as some weird code,
      // you may have to turn on Unicode encoding for text inside your code editor. - JL
      'Choose the shape/letter on the <span style="font-weight:bold">left</span> by pressing the <span style="font-weight:bold">F</span> key <br>'+
      'or the shape/letter on the <span style="font-weight:bold">right</span> by pressing the <span style="font-weight:bold">J</span> key '+
      'on your keyboard.<br><br>',
      '<p style = "text-align: center; font-size: 28px; font-weight: bold">Press spacebar to continue.</p>']
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,200,100),
      minimum_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1000,1000,1000,500]};
        }
      }
    };
    break;

  // instructions page about the number of "games" and that certain games have multiple rounds
  case 'rounds':
  //debug : need to change the names of the different conditions
    instr_gen1 = { stimuli: [
      '<p style = "font-size: 26px; font-weight: bold; text-decoration:underline">About the experiment</p><br>', // 1
      'In this experiment, you will play on 3 slot machines,<br>', // 2
      'and you will play each slot machine twice.<br><br>', // 3
      'All three slot machines will ask you to choose between two shapes.<br>', // 4
      'After each choice, you will obtain the number of points associated with the chosen shape, from 1 to 99.<br><br>', // 5
      'Slot machines 1 and 2 will sometimes change the shapes between which you have to choose.<br>', // 6
      'Slot machine 3 will always use the same pair of shapes.', // 7
      '<p style = "font-style:italic"><span style = "font-weight:bold">Note:</span> there is no carryover from one round to the next, that is,<br>'+
      'there is nothing about any game/round that will help you in understanding the symbols of any other game/round.</p><br><br>', // 8
      'This may be a lot of information, but you will be able to practice later with examples.<br><br>', // 9
      '<p style = "text-align: center; font-size: 28px; font-weight: bold">Press spacebar to continue.</p>'] // 10
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [500,1000,1000,1000,1000,1000,1000,1000,1000,500]};
        }
      }
    };
    break;

  // instructions page about the how to use the points given by the shapes, and what it means to win
  case 'points1':
    instr_gen1 = { stimuli: [
      "One of the shapes (the ‘good’ shape) gives on average<br>",
      'more points than the other (the ‘bad’ shape),<br>',
      "but you won’t be told which is which:<br>",
      "you will have to find out by playing the slot machine.<br><br>",
      "What makes your task challenging is that<br>",
      "the bad shape can sometimes give more points than the good shape, and vice versa.<br><br>",
      'This means that you can’t fully trust a single choice to know which shape is the good one.<br><br>',
      '<span style = "font-weight: bold">This is important because your performance will be assessed not by the number of points you have received,<br>',
      '<span style = "font-weight: bold">but by the number of times you have chosen the good shape.</p><br><br>',
      '<span style = "font-weight: bold; text-style: italic">Note: You will be awarded a bonus for good performance on the task.</span><br><br>',
      '<p style = "font-size: 24px; font-weight: bold">Press spacebar to continue.</p>']
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,500]};
        }
      }
    };
    break;

  case 'points2':
    instr_gen1 = { stimuli: [
      'The good shape gives more than 50 points on average,<br>',
      'whereas the bad shape gives less than 50 points on average.<br><br>',
      'It will be your job to find out which is which,<br>',
      'and to choose the good shape as much as possible.<br><br>',
      'You will now play <span style = "font-weight:bold">slot machine 1</span> to get some training with the task.<br><br>',
      'Remember that this slot machine will sometimes change the shapes between which you have to choose. <br><br>',
      'When this happens, you will have to find out which of the two new shapes is the good shape.<br><br>',
      '<p style = "font-size: 24px; font-weight: bold">Press spacebar to continue.</p>']
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1000,1000,1000,1000,1000,1000,1000,500]};
        }
      }
    };
    break;

  // REF prompt needs to present the idea that one should not be overly influenced by the fluctuations once one is sure of the source
  case 'ref_instructions':
    instr_gen1 = { stimuli: [
      'Before we begin the first game, let\'s do some training.<br><br>',
      'You will now play a short version of <span style = "font-weight:bold">slot machine 1</span>.<br><br>',
      'Are you ready?<br><br>',
      '<p style = "text-align: center; font-size: 28px; font-weight: bold">Press spacebar to continue.</p>']
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1500,1000,1000,500]};
        }
      }
    };
    break;

  // VOL prompt needs to present the idea to keep in mind that the fluctuations are indeed informative
  // also, are participants clued in on the continuous nature of the fluctuations?
  case 'vol_instructions':
    instr_gen1 = { stimuli: [
      'You will now play <span style = "font-weight:bold">slot machine 3</span>,<br>',  // 1
      'which differs from slot machines 1 and 2<br>', // 2
      'because it always uses the same pair of shapes which are colored circles<br>', // 3
      'rather than the symbols used by slot machines 1 and 2.<br><br>', // 4
      'As in slot machines 1 and 2, on each trial there is always one good shape and one bad shape.<br><br>', // 5
      'But the difference is that the better shape and the bad shape<br>', // 6
      'will sometimes switch from one to the other at random points during the game.<br><br>', // 7
      'This means that you will need to keep looking out for these switches,<br>', // 8
      'because your performance will be assessed by the number of times you have chosen the good shape.<br><br>', // 9
      '<span style = "text-style: italic"><span style ="font-weight: bold">Note:</span> '+
      'these switches will not occur one immediately after another, but there will be multiple switches during each game.</span><br><br>', // 10
      'As before, let\'s first train a on shorter version of <span style = "font-weight:bold">Slot machine 3</span>.<br><br>', // 11
      '<p style = "text-align: center; font-size: 28px; font-weight: bold">Press spacebar to continue.</p>'] // 12
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,500]};
        }
      }
    };
    break;

  case 'unp_instructions':
    instr_gen1 = { stimuli: ['You will now play <span style = "font-weight:bold">slot machine 2</span>,<br>',
      'which works exactly like slot machine 1 except<br>',
      'that the difference between the number of points associated<br>',
      'with the good and bad shapes is less obvious.<br><br>',
      'It will thus be more challenging to find out which shape is the good shape.<br><br>',
      '<p style = "text-align: center; font-size: 28px; font-weight: bold">Press spacebar to continue.</p>']
    };
    instructions_var = {
      type: 'html-keyboard-response-sequential-faded',
      stimulus: instr_gen1,
      choices: [' '],
      fadein_duration: duration_fn(debugFlag,1000,100),
      fadeout_duration: duration_fn(debugFlag,1000,100),
      individual_durations: () => {
        if (debugFlag) {
            return {durations: []};
        } else {
            return {durations: [1500,1500,1500,1500,1500,500]};
        }
      }
    };
    break;

}

return instructions_var;
}
