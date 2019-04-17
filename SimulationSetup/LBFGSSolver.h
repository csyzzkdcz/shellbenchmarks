
//bool SimulationSetupDynamicSolver::lineSearch(std::shared_ptr<SensitiveAnalysis> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
//{
//    double c1 = 0.1;
//    double c2 = 0.9;
//    double alpha = M_MIN;
//    double beta = M_MAX;
//
//    Eigen::VectorXd grad;
//    double orig = op->value(L, S, pos);
//    op->gradient(L, S, pos, grad);
//    double deriv = dir.dot(grad);
//
//    while (true)
//    {
//        Eigen::VectorXd newdE;
//        Eigen::VectorXd newL = L + rate*dir.segment(0, L.size());
//        Eigen::MatrixXd newPos = pos;
//        Eigen::VectorXd newS = S + rate*dir.segment(L.size(), S.size());
//        op->projectBack(newL, newS, newPos);
//        double newenergy = op->value(newL, newS, newPos);
//        op->gradient(newL, newS, newPos, newdE);
//
//        std::cout << "Trying rate = " << rate << ", energy now " << newenergy<<", L update "<<(newL-L).norm()<<", S update "<<(newS - S).norm() << ", L grad "<<newdE.segment(0, L.rows()).norm()<<" , S grad "<<newdE.segment(L.rows(), S.rows()).norm() << std::endl;
//
//        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
//        {
//            //            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
//            beta = rate;
//            rate = 0.5*(alpha + beta);
//            if (beta - alpha < 1e-15)
//            {
//                rate = 1e-15;
//                std::cout<<"Line Search failed, finished with Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                //                pos = newPos;
//                return false;
//            }
//        }
//        else if (newdE.dot(dir) < c2*deriv)
//        {
//            //            std::cout<<"Voilate the second Wolfe Condition"<<std::endl;
//            alpha = rate;
//            if (beta == M_MAX)
//            {
//                rate = 2 * alpha;
//            }
//            else
//            {
//                rate = 0.5*(alpha + beta);
//            }
//
//            if (beta - alpha < 1e-10)
//            {
//                if(newenergy > orig)
//                {
//                    std::cout<<"Line Search failed with beta - alph < 1e-10 without decreasing energy, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//                else
//                {
//                    std::cout<<"Line Search succeed with beta - alph < 1e-10 without sufficient decrease, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//            }
//        }
//        else
//        {
//            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
//            return true;
//        }
//    }
//}


//    //////////////////////////////////////// Optimization ///////////////////////////////////////////////
//    int m = 10;
//    int DIM = L.rows() + S.rows();
//    Eigen::VectorXd grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM), z(DIM), x(DIM);
//    Eigen::MatrixXd sVector = Eigen::MatrixXd ::Zero(DIM, m);
//    Eigen::MatrixXd yVector = Eigen::MatrixXd::Zero(DIM, m);
//    Eigen::Matrix<double, Eigen::Dynamic, 1> alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m);
//    op->gradient(L, S, pos, grad);
//    Eigen::VectorXd L_old = L;
//    Eigen::VectorXd S_old = S;
//    grad_old = grad;
//    double fmin = op->value(L, S, pos);
//
//    std::cout<<"Initial Shape Difference: "<<op->computeDifference(pos)<<std::endl;
//    std::cout<<"Initial Abar Penalty: "<<op->computeAbarSmoothness(L)<<" Penalty Coefficient: "<<abarCoef<<std::endl;
//    if(selectedDynamicType == "ABbarPos")
//    {
//        auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
//        double smoothnessBbarValue = op1->computeBbarSmoothness(L, S, pos);
//        std::cout<<"Initial Bbar Penalty: "<<smoothnessBbarValue<<" Penalty Coefficient: "<<bbarCoef<<std::endl;
//    }
//
//    int iter = 0;
//    double gradNorm = 0;
//    std::cout<<"Simulation start (L-BFGS)!! Initial function value is: "<<fmin<<std::endl<<std::endl;
//    while(true)
//    {
//
//        //        const double relativeEpsilon = static_cast<double>(1e-6) * std::max<double>(static_cast<double>(1.0), L.norm());
//        //
//        //        if (grad.norm() < relativeEpsilon)
//        //            break;
//        //
//        /////////////////////////////////////////////// L-BFGS /////////////////////////////////////////////////////////////////
//        double H0k = 1;
//        q = -grad;
//        const int k = std::min<double>(m, iter);
//        // for i = k − 1, k − 2, . . . , k − m§
//        for (int i = k - 1; i >= 0; i--)
//        {
//            if(k < m)
//                break;
//            // alpha_i <- rho_i*s_i^T*q
//            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
//            //            std::cout<<rho<<std::endl;
//            alpha(i) = rho * sVector.col(i).dot(q);
//            // q <- q - alpha_i*y_i
//            q = q - alpha(i) * yVector.col(i);
//        }
//        // z <- H_k^0*q
//        // update the scaling factor
//        if(k >= m)
//        {
//            H0k = yVector.col(m-1).dot(sVector.col(m-1)) / yVector.col(m-1).dot(yVector.col(m-1));
//        }
//        std::cout<<"H0K: "<<H0k<<std::endl;
//        z = H0k * q;
//        //for i k − m, k − m + 1, . . . , k − 1
//        for (int i = 0; i <= k-1; i++)
//        {
//            if(k < m)
//                break;
//            // beta <- rho_i * y_i^T * r
//            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
//            double beta = rho * yVector.col(i).dot(z);
//            // z <- z + s_i * ( alpha_i - beta)
//            z = z + sVector.col(i) * (alpha(i) - beta);
//        }
//        std::cout<<"z: "<<z.norm()<<", z_l: "<<z.segment(0, L.size()).norm()<<", z_s: "<<z.segment(L.size(), S.size()).norm()<<std::endl;
//        x.segment(0, L.rows()) = L;
//        x.segment(L.rows(), S.rows()) = S;
//        double rate = 0;
//        if(k<m)
//        {
//            rate = std::min(x.norm() / sqrt(x.size()) * 1e-3, 1e-5/grad.norm() * x.norm());
//            std::cout<<"initial rate: "<<rate<<" , (L,S) norm : "<<x.norm()<<" , direction norm: "<<grad.norm()<<", L direction norm: "<<grad.segment(0, L.size()).norm()<<", S direction norm: "<<grad.segment(L.size(), S.size()).norm()<<std::endl;
//            //            std::cout<<x.norm() / sqrt(x.size()) * 1e-3<<" "<<1e-5/grad.norm() * x.norm()<<std::endl;
//            bool isSuccess = lineSearch(op, L, S, pos, -grad, rate);
//            L = L - rate * grad.segment(0,L.size());
//            S = S - rate * grad.segment(L.size(), S.size());
//            op->projectBack(L, S, pos);
//            op->gradient(L, S, pos, grad);
//
//            y = grad - grad_old;
//            s.segment(0, L.size()) = L - L_old;
//            s.segment(L.size(), S.size()) = S - S_old;
//            grad_old = grad;
//            L_old = L;
//            S_old = S;
//        }
//        else
//        {
//            rate = 1.0;
//            std::cout<<"initial rate: "<<rate<<" , (L,S) norm : "<<x.norm()<<" , direction norm: "<<z.norm()<<", L direction norm: "<<z.segment(0, L.size()).norm()<<", S direction norm: "<<z.segment(L.size(), S.size()).norm()<<std::endl;
//            bool isSuccess = lineSearch(op, L, S, pos, z, rate);
//            if(isSuccess)
//            {
//                std::cout<<"L-BFGS succeeded!"<<std::endl;
//                L = L + rate * z.segment(0,L.size());
//                S = S + rate * z.segment(L.size(), S.size());
//                op->projectBack(L, S, pos);
//                op->gradient(L, S, pos, grad);
//
//                y = grad - grad_old;
//                s.segment(0, L.size()) = L - L_old;
//                s.segment(L.size(), S.size()) = S - S_old;
//                grad_old = grad;
//                L_old = L;
//                S_old = S;
//            }
//            else
//            {
//                std::cout<<"L-BFGS failed, using GD instead!!"<<std::endl;
//                rate = std::min(x.norm() / sqrt(x.size()) * 1e-3, 1e-5/grad.norm() * x.norm());
//                bool isSuccess = lineSearch(op, L, S, pos, -grad, rate);
//                L = L - rate * grad.segment(0,L.size());
//                S = S - rate * grad.segment(L.size(), S.size());
//                op->projectBack(L, S, pos);
//                op->gradient(L, S, pos, grad);
//
//                y = grad - grad_old;
//                s.segment(0, L.size()) = L - L_old;
//                s.segment(L.size(), S.size()) = S - S_old;
//                grad_old = grad;
//                L_old = L;
//                S_old = S;
//            }
//        }
//        // update the history
//        if (iter < m)
//        {
//            sVector.col(iter) = s;
//            yVector.col(iter) = y;
//        }
//        else
//        {
//            sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
//            sVector.rightCols(1) = s;
//            yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
//            yVector.rightCols(1) = y;
//        }
//        std::cout<<"s.dot(y) = "<<s.dot(y)<<std::endl;
//
//        ///////////////////////////////////////////// G-N /////////////////////////////////////////////////////////////////
//        //        double rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
//        //        Eigen::SparseMatrix<double> H = (grad * grad.transpose()).sparseView();
//        //        Eigen::SparseMatrix<double> I(H.rows(), H.cols());
//        //        I.setIdentity();
//        //        H += 1e-6 * I;
//        //        Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
//        //        solver.compute(H);
//        //        Eigen::VectorXd dir = solver.solve(-grad);
//        //        bool isSuccess = lineSearch(op, L, pos, dir, rate);
//        //
//        //        L = L + rate * dir;
//        //        grad_old = grad;
//        //        projectBackOp(sff, L, pos);
//        //        op.gradient(L, pos, grad);
//        //
//
//        /////////////////////////////////////////////// GD /////////////////////////////////////////////////////////////////
//        //        // find steplength
//        //        x.segment(0, L.rows()) = L;
//        //        x.segment(L.rows(), S.rows()) = S;
//        //        double rate = std::min(x.norm() / sqrt(x.size()) * 1e-3, 1e-5/grad.norm() * x.norm());
//        //        bool isSuccess = lineSearch(op, L, S, pos, -grad, rate);
//        //
//        //        L = L - rate * grad.segment(0,L.size());
//        //        S = S - rate * grad.segment(L.size(), S.size());
//        //        op->projectBack(L, S, pos);
//        //        op->gradient(L, S, pos, grad);
//        //
//        //        y = grad - grad_old;
//        //        s.segment(0, L.size()) = L - L_old;
//        //        s.segment(L.size(), S.size()) = S - S_old;
//        //        grad_old = grad;
//        //        L_old = L;
//        //        S_old = S;
//
//        std::cout<<std::endl<< "iter: "<<iter<< ", Rate: "<<rate<< ", f = " <<  op->value(L, S, pos)<<", ||g||_inf "<< grad.lpNorm<Eigen::Infinity>()<<std::endl;
//        std::cout <<"Shape Difference = "<<op->computeDifference(pos)<<", Abar Penalty: "<<op->computeAbarSmoothness(L);
//
//        if(selectedDynamicType == "ABbarPos")
//        {
//            auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
//            double smoothnessBbarValue = op1->computeBbarSmoothness(L, S, pos);
//            std::cout<<", Bbar Penalty: "<<smoothnessBbarValue<<std::endl;
//        }
//        else
//        {
//            std::cout<<std::endl;
//        }
//
//        std::cout<<"||g_L||_inf "<<grad.segment(0, L.rows()).lpNorm<Eigen::Infinity>()<<", ||g_s||_inf "<<grad.segment(L.rows(), S.rows()).lpNorm<Eigen::Infinity>()<<std::endl;
//        std::cout<<"abar change "<<s.segment(0, L.size()).lpNorm<Eigen::Infinity>()<<", bbar changes: "<<s.segment(L.size(), S.size()).lpNorm<Eigen::Infinity>()<<", gradient change "<<y.lpNorm<Eigen::Infinity>()<<std::endl<<std::endl;
//
//
//
//        iter++;
//        gradNorm = grad.template lpNorm<Eigen::Infinity>();
//        if(isnan(gradNorm) || isnan(s.dot(y)))
//        {
//            std::cout<<"Something wrong happened!!"<<std::endl;
//            std::cout<<s.norm()<<std::endl;
//            std::cout<<y.norm()<<std::endl;
//            std::ofstream outfile("error.dat", std::ios::trunc);
//            int nverts = targetPos.rows();
//            int nfaces = mesh.nFaces();
//
//            outfile<<thickness<<"\n";
//            outfile<<abarCoef<<"\n";
//            if(selectedDynamicType != "AbarPos")
//                outfile<<bbarCoef<<"\n";
//            outfile<<smoothCoef<<"\n";
//            outfile<<3*nverts<<"\n";
//            outfile<<3*nfaces<<"\n";
//
//            std::cout<<3*nverts + 3*nfaces<<std::endl;
//
//            for(int i=0;i<nverts;i++)
//            {
//                outfile<<std::setprecision(16)<<pos(i, 0)<<"\n";
//                outfile<<std::setprecision(16)<<pos(i, 1)<<"\n";
//                outfile<<std::setprecision(16)<<pos(i, 2)<<"\n";
//            }
//
//            for(int i=0;i<3*nfaces;i++)
//            {
//                outfile<<std::setprecision(16)<<L(i)<<"\n";
//            }
//            outfile<<std::setprecision(16)<<L(L.size()-1);
//            outfile.close();
//            igl::writeOBJ("resampled_error.obj", pos, mesh.faces());
//        }
//        if(iter == MAX_ITR)
//        {
//            std::cout<<"Maximun iteration reached!"<<std::endl;
//            break;
//        }
//        if(gradNorm <= M_TOL * std::max<double>(static_cast<double>(1.0), L.norm()))
//        {
//            std::cout<<"Force norm is less than "<<M_TOL<<std::endl;
//            break;
//        }
//        if(s.template lpNorm<Eigen::Infinity>() <= M_TOL)
//        {
//            std::cout<<"variable update is less than "<<M_TOL<<std::endl;
//            break;
//        }
//        if(y.template lpNorm<Eigen::Infinity>() <= M_TOL)
//        {
//            std::cout<<"gradient update is less than "<<M_TOL<<std::endl;
//            break;
//        }
//        if(iter % 10 == 0)
//        {
//            double f = op->value(L, S, pos);
//            if(f < fmin)
//            {
//                saveAbars(L, S, pos, false);
//                fmin = f;
//            }
//        }
//    }
//
//    double f = op->value(L, S, pos);
//    if(f < fmin)
//    {
//        saveAbars(L, S, pos, false);
//        fmin = f;
//    }
