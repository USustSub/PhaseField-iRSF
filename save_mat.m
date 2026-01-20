function save_mat ( file_ , ff , inc ) 

    filename = sprintf(['out/' (file_) '/jj%02d.mat'], inc);
    save(filename, 'ff');

end