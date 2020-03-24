function labels=DIR_Process_labels(labels_search,level)
labels.sx=labels_search(level,1);
labels.sy=labels_search(level,2);
labels.sz=labels_search(level,3);
labels.n_each_layer=labels.sx*labels.sy;
labels.nlabels=labels.sx*labels.sy*labels.sz;
labels.hz=floor(labels.sz/2);
labels.hx=floor(labels.sx/2);
labels.hy=floor(labels.sy/2);

end