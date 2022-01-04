function [Bi] = Bidiag_Francis_Step(Bi)    
    tau10 = Bi(1, 2) * Bi(1,1);
    tau00 = Bi(1, 1)^2;
    tau_last = Bi(end-1,end)^2 + Bi(end, end)^2;
    G0 = Givens_rotation([tau00 - tau_last; tau10]);
    Bi(1:2,1:2) = Bi(1:2,1:2) * G0;
        
    op_on_rows = true;
    bulge_row = 2;
    bulge_col = 1;
    m = size(Bi, 1);
    n_iters = 2 * (m-1) - 1;
    for i = 1:n_iters
        if op_on_rows
            prev_row = bulge_row - 1;
            if i == n_iters
                nxt_col = bulge_col + 1;
            else
                nxt_col = bulge_col + 2;
            end
            G = Givens_rotation(Bi(prev_row:bulge_row, bulge_col));
            Bi(prev_row:bulge_row,bulge_col:nxt_col) = G' * ...
                Bi(prev_row:bulge_row,bulge_col:nxt_col);
            bulge_row = prev_row;
            bulge_col = nxt_col;
        else
            prev_col = bulge_col - 1;
            nxt_row = bulge_row + 2;
            G = Givens_rotation(Bi(bulge_row, prev_col:bulge_col)');
            Bi(bulge_row:nxt_row, prev_col:bulge_col) = ...
                Bi(bulge_row:nxt_row, prev_col:bulge_col) * G;
            bulge_row = nxt_row;
            bulge_col = prev_col;
        end
        op_on_rows = ~op_on_rows;
    end
end