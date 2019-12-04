import mido
print mido.get_input_names()
#msg= mido.Message('note_on', note=60)

last_midi_event={}
inport= mido.open_input('Network Session 1')
#inport= mido.open_input('Launchkey MK2 25 Launchkey MIDI')

# non blocking input
while True:
    for msg in inport.iter_pending():
        print msg
        # assuming we get only cc changes
        tmp = "%s" % (msg)
        tmp=tmp.split(' ')
        # channel third
        chan=tmp[2].split('=')[1]
        print chan
        val=tmp[3].split('=')[1]
        last_midi_event[chan]=val

