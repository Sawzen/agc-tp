def get_chunks(sequence, chunk_size):
	list_seg =[]


	for i in (range(0, len(sequence) , chunk_size)):
		if i+chunk_size<=len(sequence):
			list_seg.append(sequence[i:i+chunk_size])
	if len(list_seg)>=4:
		return list_seg


