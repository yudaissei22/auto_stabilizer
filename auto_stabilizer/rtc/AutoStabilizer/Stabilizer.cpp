#include "Stabilizer.h"
#include "MathUtil.h"
#include <cnoid/Jacobian>
#include <cnoid/EigenUtil>
#include <cnoid/src/Body/InverseDynamics.h>

void Stabilizer::initStabilizerOutput(const GaitParam& gaitParam,
                                      cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Vector3& o_stTargetZmp, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage) const{

  o_stOffsetRootRpy.reset(cnoid::Vector3::Zero());
  o_stTargetZmp = gaitParam.refZmpTraj[0].getStart();
  for(int i=0;i<o_stServoPGainPercentage.size();i++){
    o_stServoPGainPercentage[i].reset(100.0);
    o_stServoDGainPercentage[i].reset(100.0);
  }
}

bool Stabilizer::execStabilizer(const GaitParam& gaitParam, double dt, bool useActState,
                                cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose, cnoid::Vector3& o_stTargetZmp, std::vector<cnoid::Vector6>& o_stEETargetWrench, cnoid::Vector3& o_stTargetCogAcc) const{
  // - root attitude control
  // - 現在のactual重心位置から、目標ZMPを計算
  // - 目標ZMPを満たすように目標足裏反力を計算
  // - 目標反力を満たすように重力補償+仮想仕事の原理

  // masterのhrpsysではSwingEEModificationを行っていたが、地面についているときに、time_constでもとの位置に戻るまでの間、足と重心の相対位置が着地位置タイミング修正計算で用いるものとずれることによる着地位置修正パフォーマンスの低下のデメリットの方が大きいので、削除した
  // masterのhrpsysやもとのauto_stabilizerでは位置制御+DampingControlをサポートしていたが、位置制御+DampingControlは実機での目標接触力への追従に遅れがある. FootGuidedControlでは、目標ZMPの位相を進めたり、ZMPの追従遅れを考慮した目標ZMP計算を行ったりをしていないので、遅れに弱い. そのため、位置制御+DampingControlは削除し、所謂TorqueControlのみをサポートしている.

  // root attitude control
  this->moveBasePosRotForBodyRPYControl(dt, gaitParam, useActState,// input
                                        o_stOffsetRootRpy, o_stTargetRootPose); // output

  // 現在のactual重心位置から、目標ZMPを計算
  cnoid::Vector3 tgtForce; // generate frame
  this->calcZMP(gaitParam, dt, useActState, // input
                o_stTargetZmp, tgtForce, o_stTargetCogAcc); // output

  // 目標ZMPを満たすように目標EndEffector反力を計算
  this->calcWrench(gaitParam, o_stTargetZmp, tgtForce, useActState,// input
                   o_stEETargetWrench); // output

  return true;
}

bool Stabilizer::moveBasePosRotForBodyRPYControl(double dt, const GaitParam& gaitParam, bool useActState,
                                                 cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const{

  // stOffsetRootRpyを計算
  if(!useActState){
      o_stOffsetRootRpy.interpolate(dt);
  }else{
    cnoid::Vector3 stOffsetRootRpy = o_stOffsetRootRpy.value(); // gaitParam.footMidCoords frame

    cnoid::Matrix3 rootRErrorGenerateFrame = gaitParam.refRobot->rootLink()->R() * gaitParam.actRobot->rootLink()->R().transpose(); // generate frame
    cnoid::Matrix3 rootRError = gaitParam.footMidCoords.value().linear().transpose() * rootRErrorGenerateFrame/*generate frame*/ * gaitParam.footMidCoords.value().linear(); // gaitParam.footMidCoords frame
    cnoid::Vector3 rootRpyError = cnoid::rpyFromRot(rootRError); // gaitParam.footMidCoords frame

    for (size_t i = 0; i < 2; i++) {
      stOffsetRootRpy[i] += (this->bodyAttitudeControlGain[i] * rootRpyError[i] - 1.0/this->bodyAttitudeControlTimeConst[i] * stOffsetRootRpy[i]) * dt;
      stOffsetRootRpy[i] = mathutil::clamp(stOffsetRootRpy[i], this->bodyAttitudeControlCompensationLimit[i]);
    }
    stOffsetRootRpy[2] = 0.0;

    o_stOffsetRootRpy.reset(stOffsetRootRpy);
  }

  // stTargetRootPoseを計算
  o_stTargetRootPose.translation() = gaitParam.refRobot->rootLink()->p();
  o_stTargetRootPose.linear() /*generate frame*/= gaitParam.footMidCoords.value().linear() * cnoid::rotFromRpy(o_stOffsetRootRpy.value()/*gaitParam.footMidCoords frame*/) * gaitParam.footMidCoords.value().linear().transpose() * gaitParam.refRobot->rootLink()->R()/*generate frame*/;
  return true;
}

bool Stabilizer::calcZMP(const GaitParam& gaitParam, double dt, bool useActState,
                         cnoid::Vector3& o_tgtZmp, cnoid::Vector3& o_tgtForce, cnoid::Vector3& o_tgtCogAcc) const{
  cnoid::Vector3 cog = useActState ? gaitParam.actCog : gaitParam.genCog;
  cnoid::Vector3 cogVel = useActState ? gaitParam.actCogVel.value() : gaitParam.genCogVel;
  cnoid::Vector3 DCM = cog + cogVel / gaitParam.omega;
  const std::vector<cnoid::Position>& EEPose = useActState ? gaitParam.actEEPose : gaitParam.abcEETargetPose;

  cnoid::Vector3 tgtZmp;
  if(gaitParam.footstepNodesList[0].isSupportPhase[RLEG] || gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    tgtZmp = footguidedcontroller::calcFootGuidedControl(gaitParam.omega,gaitParam.l,DCM,gaitParam.refZmpTraj);
    if(tgtZmp[2] >= gaitParam.actCog[2]) tgtZmp = gaitParam.actCog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0); // 下向きの力は受けられないので
    else{
      // truncate zmp inside polygon. actual robotの関節角度を用いて計算する
      std::vector<cnoid::Vector3> vertices; // generate frame. 支持点の集合
      for(int i=0;i<NUM_LEGS;i++){
        if(!gaitParam.footstepNodesList[0].isSupportPhase[i]) continue;
        for(int j=0;j<gaitParam.legHull[i].size();j++){
          cnoid::Vector3 p = EEPose[i]*gaitParam.legHull[i][j];
          if(p[2] > gaitParam.actCog[2] - 1e-2) p[2] = gaitParam.actCog[2] - 1e-2; // 重心よりも支持点が高いと射影が破綻するので 
          vertices.push_back(p);
        }
      }
      tgtZmp = mathutil::calcInsidePointOfPolygon3D(tgtZmp,vertices,gaitParam.actCog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0));
      // TODO. 角運動量オフセット.
    }
  }else{ // 跳躍期
    tgtZmp = cog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0);
  }
  cnoid::Vector3 tgtCog,tgtCogVel,tgtCogAcc,tgtForce;
  footguidedcontroller::updateState(gaitParam.omega,gaitParam.l,cog,cogVel,tgtZmp,gaitParam.genRobot->mass(),dt,
                                      tgtCog, tgtCogVel, tgtCogAcc, tgtForce);

  // tgtForceにrefEEWrenchのXY成分を足す TODO

  o_tgtZmp = tgtZmp;
  o_tgtForce = tgtForce;
  o_tgtCogAcc = tgtCogAcc;
  return true;
}

bool Stabilizer::calcWrench(const GaitParam& gaitParam, const cnoid::Vector3& tgtZmp/*generate座標系*/, const cnoid::Vector3& tgtForce/*generate座標系 ロボットが受ける力*/, bool useActState,
                            std::vector<cnoid::Vector6>& o_tgtEEWrench) const{
  std::vector<cnoid::Vector6> tgtEEWrench(gaitParam.eeName.size(), cnoid::Vector6::Zero()); /* 要素数EndEffector数. generate frame. EndEffector origin*/

  for(int i = 0;i<gaitParam.eeName.size();i++){
    tgtEEWrench[i] = gaitParam.refEEWrench[i];
  }

  /*
    legは、legから受けるwrenchの和がtgtZmp, tgtForceを満たすように.
    非Support期のlegには分配せずゼロを入れる. 全てのlegが非Support期なら分配計算すら行わない
    actual robotの関節角度を用いて計算する
    各EEFのwrenchを、polygonの各頂点からのSPAN表現で考える.
    各頂点のfx, fy, fzの向きは、合力の向きと同じで、ノルムだけを変数とする. 合力がtgtForce, ZMPがtgtZmpになるように、ノルムの値を求める.
      - nzが反映できない、力の向きの冗長性を利用できない、摩擦係数を考慮できない、といった欠点がある. 二次元動歩行なので、まずは物理的・数学的厳密性や冗長性の利用よりもシンプルさ、ロバストさを優先する. そのあたりをこだわりたいなら三次元多点接触でやる.
      - 動歩行の途中の一歩で偶然actualの足が90度以上倒れた姿勢で地面につくことがあるので、そうなったときにても破綻しないことが重要.
    最後に、FACE表現に変換する.

    階層QPのタスクは次の通り
    変数: SPAN表現のノルム. 0~1
    1. ノルム>0. 合力がtgtForce.
    2. ZMPがtgtZmp
    3. 各脚の各頂点のノルムの重心がCOPOffsetと一致 (fzの値でスケールされてしまうので、alphaを用いて左右をそろえる)
    4. ノルムの2乗和の最小化 (3の中で微小な重みで一緒にやる)
  */
  // 計算時間は、tgtZmpが支持領域内に無いと遅くなるなので、事前に支持領域内に入るように修正しておくこと
  const std::vector<cnoid::Position>& EEPose = useActState ? gaitParam.actEEPose : gaitParam.abcEETargetPose;

  if(gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && !gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    if(gaitParam.isManualControlMode[LLEG].getGoal() == 0.0) tgtEEWrench[LLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
    tgtEEWrench[RLEG].head<3>() = tgtForce;
    tgtEEWrench[RLEG].tail<3>() = (tgtZmp - EEPose[RLEG].translation()).cross(tgtForce);
  }else if(!gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    if(gaitParam.isManualControlMode[RLEG].getGoal() == 0.0) tgtEEWrench[RLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
    tgtEEWrench[LLEG].head<3>() = tgtForce;
    tgtEEWrench[LLEG].tail<3>() = (tgtZmp - EEPose[LLEG].translation()).cross(tgtForce);
  }else if(!gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && !gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    if(gaitParam.isManualControlMode[RLEG].getGoal() == 0.0) tgtEEWrench[RLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
    if(gaitParam.isManualControlMode[LLEG].getGoal() == 0.0) tgtEEWrench[LLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
  }else if(tgtForce.norm() == 0){
    if(gaitParam.isManualControlMode[RLEG].getGoal() == 0.0) tgtEEWrench[RLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
    if(gaitParam.isManualControlMode[LLEG].getGoal() == 0.0) tgtEEWrench[LLEG].setZero(); // Manual Control ModeであればrefEEWrenchをそのまま使う
  }else{
    int dim = gaitParam.legHull[RLEG].size() + gaitParam.legHull[LLEG].size();
    {
      // 1. ノルム>0. 合力がtgtForce.

      // 合力がtgtForce. (合計が1)
      this->constraintTask_->A() = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,dim);
      for(int i=0;i<dim;i++) this->constraintTask_->A().insert(0,i) = 1.0;
      this->constraintTask_->b() = Eigen::VectorXd::Ones(1);
      this->constraintTask_->wa() = cnoid::VectorX::Ones(1);

      // 各値が0~1
      this->constraintTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(dim,dim);
      for(int i=0;i<dim;i++) this->constraintTask_->C().insert(i,i) = 1.0;
      this->constraintTask_->dl() = Eigen::VectorXd::Zero(dim);
      this->constraintTask_->du() = Eigen::VectorXd::Ones(dim);
      this->constraintTask_->wc() = cnoid::VectorX::Ones(dim);

      this->constraintTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->constraintTask_->toSolve() = false;
      this->constraintTask_->settings().verbose = 0;
    }
    {
      // 2. ZMPがtgtZmp
      // tgtZmpまわりのトルクの和を求めて、tgtForceの向きの単位ベクトルとの外積が0なら良い
      //this->tgtZmpTask_->A() = Eigen::SparseMatrix<double,Eigen::RowMajor>(3,dim);
      cnoid::Vector3 tgtForceDir = tgtForce.normalized();
      int idx = 0;
      Eigen::SparseMatrix<double,Eigen::ColMajor> A_ColMajor(3,dim); // insert()する順序がColMajorなので、RowMajorのAに直接insertすると計算効率が著しく悪い(ミリ秒単位で時間がかかる).
      for(int i=0;i<NUM_LEGS;i++){
        for(int j=0;j<gaitParam.legHull[i].size();j++){
          cnoid::Vector3 pos = EEPose[i] * gaitParam.legHull[i][j];
          cnoid::Vector3 a = tgtForceDir.cross( (pos - tgtZmp).cross(tgtForce));
          for(int k=0;k<3;k++) A_ColMajor.insert(k,idx) = a[k];
          idx ++;
        }
      }
      this->tgtZmpTask_->A() = A_ColMajor;
      this->tgtZmpTask_->b() = Eigen::VectorXd::Zero(3);
      this->tgtZmpTask_->wa() = cnoid::VectorX::Ones(3);

      this->tgtZmpTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->tgtZmpTask_->dl() = Eigen::VectorXd::Zero(0);
      this->tgtZmpTask_->du() = Eigen::VectorXd::Ones(0);
      this->tgtZmpTask_->wc() = cnoid::VectorX::Ones(0);

      this->tgtZmpTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->tgtZmpTask_->toSolve() = false; // 常にtgtZmpが支持領域内にあるなら解く必要がないので高速化のためfalseにする. ない場合があるならtrueにする. calcWrenchでtgtZmpをtruncateしているのでfalseでよい
      this->tgtZmpTask_->settings().verbose = 0;
    }
    {
      // 3. 各脚の各頂点のノルムの重心がCOPOffsetと一致 (fzの値でスケールされてしまうので、alphaを用いて左右をそろえる)

      // 各EndEffectorとtgtZmpの距離を用いてalphaを求める
      std::vector<double> alpha(NUM_LEGS);
      {
        cnoid::Vector3 rleg2leg = EEPose[LLEG].translation() - EEPose[RLEG].translation();
        rleg2leg[2] = 0.0;
        if(rleg2leg.norm() == 0.0){
          alpha[RLEG] = alpha[LLEG] = 0.5;
        }else{
          cnoid::Vector3 rleg2legDir = rleg2leg.normalized();
          double rleg2llegDistance = rleg2leg.norm();
          double rleg2tgtZmpRatio = rleg2legDir.dot(tgtZmp - EEPose[RLEG].translation()) / rleg2llegDistance;
          alpha[RLEG] = mathutil::clamp(1.0 - rleg2tgtZmpRatio, 0.05, 1.0-0.05);
          alpha[LLEG] = mathutil::clamp(rleg2tgtZmpRatio, 0.05, 1.0-0.05);
        }
      }
      Eigen::SparseMatrix<double,Eigen::ColMajor> A_ColMajor(3*NUM_LEGS,dim); // insert()する順序がColMajorなので、RowMajorのAに直接insertすると計算効率が著しく悪い(ミリ秒単位で時間がかかる).
      int idx = 0;
      for(int i=0;i<NUM_LEGS;i++) {
        cnoid::Vector3 cop = EEPose[i].translation() + EEPose[i].linear() * gaitParam.copOffset[i].value();
        for(int j=0;j<gaitParam.legHull[i].size();j++){
          cnoid::Vector3 pos = EEPose[i].translation() + EEPose[i].linear() * gaitParam.legHull[i][j];
          cnoid::Vector3 a = (pos - cop) / alpha[i];
          for(int k=0;k<3;k++) A_ColMajor.insert(i*3+k,idx) = a[k];
          idx ++;
        }
      }
      this->copTask_->A() = A_ColMajor;
      this->copTask_->b() = Eigen::VectorXd::Zero(3*NUM_LEGS);
      this->copTask_->wa() = cnoid::VectorX::Ones(3*NUM_LEGS);

      this->copTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->copTask_->dl() = Eigen::VectorXd::Zero(0);
      this->copTask_->du() = Eigen::VectorXd::Ones(0);
      this->copTask_->wc() = cnoid::VectorX::Ones(0);

      this->copTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->copTask_->toSolve() = true;
      // this->copTask_->options().setToReliable();
      // this->copTask_->options().printLevel = qpOASES::PL_NONE; // PL_HIGH or PL_NONE
      this->copTask_->settings().check_termination = 5; // default 25. 高速化
      this->copTask_->settings().verbose = 0;
    }

    std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks{this->constraintTask_,this->tgtZmpTask_,this->copTask_};
    cnoid::VectorX result;
    if(!prioritized_qp_base::solve(tasks,
                                   result,
                                   0 // debuglevel
                                   )){
      // QP fail. 適当に1/2 tgtForceずつ分配してもよいが、QPがfailするのはだいたい転んでいるときなので、ゼロを入れたほうが安全
      tgtEEWrench[RLEG].setZero();
      tgtEEWrench[LLEG].setZero();
    }else{
      int idx = 0;
      for(int i=0;i<NUM_LEGS;i++){
        tgtEEWrench[i].setZero();
        for(int j=0;j<gaitParam.legHull[i].size();j++){
          tgtEEWrench[i].head<3>() += tgtForce * result[idx];
          tgtEEWrench[i].tail<3>() += (EEPose[i].linear() * gaitParam.legHull[i][j]).cross(tgtForce * result[idx]);
          idx ++;
        }
      }
    }
  }

  o_tgtEEWrench = tgtEEWrench;
  return true;
}

bool Stabilizer::calcTorque(double dt, const GaitParam& gaitParam, bool useActState, cnoid::BodyPtr& actRobotTqc, const std::vector<cnoid::Vector6>& tgtEEWrench /* 要素数EndEffector数. generate座標系. EndEffector origin*/,
                            std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage, Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd) const{

  if(!useActState){
    for(int i=0;i<actRobotTqc->numJoints();i++) actRobotTqc->joint(i)->u() = 0.0;
    for(int i=0;i<actRobotTqc->numJoints();i++){
      o_stServoPGainPercentage[i].interpolate(dt);
      o_stServoDGainPercentage[i].interpolate(dt);
    }
  } else {
    cnoid::VectorX tau_g = cnoid::VectorXd::Zero(actRobotTqc->numJoints()); // 重力
    cnoid::VectorX tau_lip = cnoid::VectorXd::Zero(actRobotTqc->numJoints()); // 倒立振子のtgtWrench
    cnoid::VectorX tau_ee = cnoid::VectorXd::Zero(actRobotTqc->numJoints()); // endEffector

    // 速度・加速度を考慮しない重力補償
    {
      actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
      actRobotTqc->rootLink()->v() = cnoid::Vector3::Zero();
      actRobotTqc->rootLink()->w() = cnoid::Vector3::Zero();
      actRobotTqc->rootLink()->dv() = cnoid::Vector3(0.0,0.0,gaitParam.g);
      actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
      for(int i=0;i<actRobotTqc->numJoints();i++){
	actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	actRobotTqc->joint(i)->dq() = 0.0;
	actRobotTqc->joint(i)->ddq() = 0.0;
      }
      actRobotTqc->calcForwardKinematics(true, true);
      cnoid::calcInverseDynamics(actRobotTqc->rootLink()); // actRobotTqc->joint()->u()に書き込まれる
      for(int i=0;i<actRobotTqc->numJoints();i++){
	tau_g[i] = actRobotTqc->joint(i)->u();
	actRobotTqc->joint(i)->u() = 0.0;
      }
    }

    // tgtEEWrench
    {
      for(int i=0;i<gaitParam.eeName.size();i++){
	cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	cnoid::MatrixXd J = cnoid::MatrixXd::Zero(6,jointPath.numJoints()); // generate frame. endeffector origin
	cnoid::setJacobian<0x3f,0,0,true>(jointPath,actRobotTqc->link(gaitParam.eeParentLink[i]),gaitParam.eeLocalT[i].translation(), // input
					  J); // output
	cnoid::VectorX tau = - J.transpose() * tgtEEWrench[i];
	for(int j=0;j<jointPath.numJoints();j++){
	  jointPath.joint(j)->u() += tau[j];
	}
      }
      for(int i=0;i<actRobotTqc->numJoints();i++){
	tau_lip[i] = actRobotTqc->joint(i)->u();
	actRobotTqc->joint(i)->u() = 0.0;
      }
    }

    {
      this->eeComTask_->A() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 * gaitParam.eeName.size() + 3,6 + actRobotTqc->numJoints());
      std::vector<Eigen::Triplet<double> > eeComTripletList_A;
      eeComTripletList_A.reserve(500);//適当
      this->eeComTask_->b() = Eigen::VectorXd::Zero(6 * gaitParam.eeName.size() + 3);
      this->eeComTask_->C() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 + actRobotTqc->numJoints(),6 + actRobotTqc->numJoints());
      std::vector<Eigen::Triplet<double> > eeComTripletList_C;
      // エンドエフェクタ分解加速度制御
      // ee_acc = ee_acc_ref + K (ee_p_ref - ee_p_act) + D (ee_vel_ref - ee_vel_act)
      // ee_vel = J * dq
      // ee_acc = J * ddq + dJ * dq
      // ここで dJ * dq はddq = 0 としたときのee_accだから、q, dqをもとにForwardKinematicsで求められる
      // J * ddq = ee_acc - dJ * dq を満たすddqを求めて、actualのq, dqと合わせてInverseDynamicsでuを求める
      
      {
	std::vector<cnoid::Vector6> ee_acc;
	ee_acc.resize(gaitParam.eeName.size());
	// ee_act
	{
	  // ee_act_ref
	  for(int i=0;i<gaitParam.eeName.size();i++){
	    ee_acc[i].head<3>() = (gaitParam.abcEETargetPose[i].translation() - 2 * eeTargetPosed[i].translation() + eeTargetPosedd[i].translation()) / dt / dt;
	    ee_acc[i].tail<3>() = (cnoid::rpyFromRot(gaitParam.abcEETargetPose[i].linear() * eeTargetPosed[i].linear().transpose()) - cnoid::rpyFromRot(eeTargetPosed[i].linear() * eeTargetPosedd[i].linear().transpose())) / dt / dt; // TODO
	  }

	  // K (ee_p_ref - ee_p_act) + D (ee_vel_ref - ee_vel_act)
	  for(int i=0;i<gaitParam.eeName.size();i++){
	    cnoid::Matrix3 eeR = gaitParam.actEEPose[i].linear();
	    cnoid::Vector6 eePoseDiffLocal; // endEfector frame
	    eePoseDiffLocal.head<3>() = eeR.transpose() * (gaitParam.abcEETargetPose[i].translation() - gaitParam.actEEPose[i].translation());
	    eePoseDiffLocal.tail<3>() = cnoid::rpyFromRot(gaitParam.actEEPose[i].linear().transpose()*gaitParam.abcEETargetPose[i].linear());
	    cnoid::Vector6 eeVelDiffLocal = (eePoseDiffLocal - eePoseDiff_prev[i]) / dt;
	    cnoid::Vector6 eePoseDiffGainLocal;
	    cnoid::Vector6 eeVelDiffGainLocal;
	    for(int j=0;j<6;j++){
	      eePoseDiffGainLocal[j] = this->ee_K[i][j] * eePoseDiffLocal[j];
	      eeVelDiffGainLocal[j] = this->ee_D[i][j] * eeVelDiffLocal[j];
	    }
	    ee_acc[i].head<3>() += eeR * eePoseDiffGainLocal.head<3>(); // generate frame
	    ee_acc[i].tail<3>() += eeR * eePoseDiffGainLocal.tail<3>(); // generate frame
	    ee_acc[i].head<3>() += eeR * eeVelDiffGainLocal.head<3>(); // generate frame
	    ee_acc[i].tail<3>() += eeR * eeVelDiffGainLocal.tail<3>(); // generate frame
	    eePoseDiff_prev[i] = eePoseDiffLocal;
	  }
	} // ee_act

	std::vector<cnoid::Vector6> dJdq;
	dJdq.resize(gaitParam.eeName.size());
	// dJ * dqを求める
	{
	  actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	  actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	  actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	  actRobotTqc->rootLink()->dv() = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
	  for(int i=0;i<actRobotTqc->numJoints();i++){
	    actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	    // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	    // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	    actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	    actRobotTqc->joint(i)->ddq() = 0.0;
	  }
	  actRobotTqc->calcForwardKinematics(true, true); // dJ * dq がactRobotTqc->link(gaitParam.eeParentLink[i])->dv(), dw()に格納される
	  actRobotTqc->calcCenterOfMass();
	  for(int i=0;i<gaitParam.eeName.size();i++){
	    dJdq[i][0] = actRobotTqc->link(gaitParam.eeParentLink[i])->dv()[0];
	    dJdq[i][1] = actRobotTqc->link(gaitParam.eeParentLink[i])->dv()[1];
	    dJdq[i][2] = actRobotTqc->link(gaitParam.eeParentLink[i])->dv()[2];
	    dJdq[i][3] = actRobotTqc->link(gaitParam.eeParentLink[i])->dw()[0];
	    dJdq[i][4] = actRobotTqc->link(gaitParam.eeParentLink[i])->dw()[1];
	    dJdq[i][5] = actRobotTqc->link(gaitParam.eeParentLink[i])->dw()[2];
	  }
	} // dJ * dq

	{ // eeComTaskを作る
	  for(int i=0;i<gaitParam.eeName.size();i++){

	    // AはJacobian
	    cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	    cnoid::MatrixXd J = cnoid::MatrixXd::Zero(6,jointPath.numJoints()); // generate frame. endeffector origin
	    cnoid::setJacobian<0x3f,0,0,true>(jointPath,actRobotTqc->link(gaitParam.eeParentLink[i]),gaitParam.eeLocalT[i].translation(), // input
					      J); // output

	    // 該当する箇所に代入
	    for (int j=0;j<jointPath.numJoints();j++) {
	      for(int k=0;k<6;k++){
		eeComTripletList_A.push_back(Eigen::Triplet<double>(k + 6*i,6+jointPath.joint(j)->jointId(),J(k,j))); // insertすると時間がかかる
	      }
	    }
	  
	    for(int j=0;j<6;j++) eeComTripletList_A.push_back(Eigen::Triplet<double>(j + 6*i,j,1.0)); // insertすると時間がかかる

	    // bはee_acc - dJ * dq
	    this->eeComTask_->b().block(i*6,0,6,1) = ee_acc[i] - dJdq[i];
	  
	  }

	  /*	  for(int i=0;i<3;i++) {
	    bool sign = gaitParam.genRobot->rootLink()->p()[i] > gaitParam.actRobot->rootLink()->p()[i]; 
	    eeComTripletList_C.push_back(Eigen::Triplet<double>(i,i,sign? 1.0 : -1.0));
	  }

	  for(int i=0;i<3;i++) {
	    bool sign = cnoid::rpyFromRot(gaitParam.genRobot->rootLink()->R())[i] > cnoid::rpyFromRot(gaitParam.actRobot->rootLink()->R())[i];
	    eeComTripletList_C.push_back(Eigen::Triplet<double>(3+i,3+i,sign? 1.0 : -1.0));
	    }*/
	
	  for(int i=0;i<actRobotTqc->numJoints();i++){
	    bool sign = gaitParam.genRobot->joint(i)->q() > gaitParam.actRobot->joint(i)->q();
	    eeComTripletList_C.push_back(Eigen::Triplet<double>(6+i,6+i,sign? 1.0 : -1.0)); // insertすると時間がかかる
	  }
	} // eeTask
      } // ee 分解加速度制御

      // 重心分解加速度制御
      // com_acc = com_acc_ref + K (com_p_ref - com_p_act) + D (com_vel_ref - com_vel_act)
      // com_vel = J * dq
      // com_acc = J * ddq + dJ * dq
      // ここで dJ * dq はddq = 0 としたときのcom_accだから、q, dqをもとにForwardKinematicsで求められる
      // J * ddq = com_acc - dJ * dq を満たすddqを求めて、actualのq, dqと合わせてInverseDynamicsでuを求める
      {
	cnoid::Vector3 com_acc;
	// com_act
	{
	  cnoid::Vector3 actCog = actRobotTqc->centerOfMass();
	  for(int i=0;i<3;i++){
	    com_acc[i] = gaitParam.stTargetCogAcc[i] + this->com_K[i] * (gaitParam.genCog[i] - actCog[i]) + this->com_D[i] * (gaitParam.genCogVel[i] - gaitParam.actCogVel.value()[i]); // TODO genCogではなくstで出されたcogをつかうこと
	  }
	}

	cnoid::Vector3 dJdq;
	// dJ * dqを求める
	{
	  dJdq = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	  actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	  actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	  actRobotTqc->rootLink()->dv() = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
	  for(int i=0;i<actRobotTqc->numJoints();i++){
	    actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	    // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	    // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	    actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	    actRobotTqc->joint(i)->ddq() = 0.0;
	  }
	  actRobotTqc->calcForwardKinematics(true, true);
	  actRobotTqc->calcCenterOfMass();
	  cnoid::Vector6 virtualWrench = cnoid::calcInverseDynamics(actRobotTqc->rootLink());
	  dJdq = virtualWrench.head<3>() / gaitParam.actRobot->mass();
	}

	{ // comTaskを作る
	  cnoid::MatrixXd CMJ = cnoid::MatrixXd::Zero(3,6 + actRobotTqc->numJoints());
	  cnoid::calcCMJacobian(actRobotTqc,nullptr,CMJ); //仕様でrootは後ろにつくので注意
	  for (int i=0;i<actRobotTqc->numJoints();i++) {
	    for(int j=0;j<3;j++) {
	      eeComTripletList_A.push_back(Eigen::Triplet<double>(6 * gaitParam.eeName.size() + j,6+i,CMJ(j,i)));
	    }
	  }
	  for (int i=0;i<6;i++) {
	    for(int j=0;j<3;j++) {
	      eeComTripletList_A.push_back(Eigen::Triplet<double>(6 * gaitParam.eeName.size() + j,i,CMJ(j,i + actRobotTqc->numJoints())));
	    }
	  }
	  this->eeComTask_->b().tail<3>() = com_acc - dJdq;
	}
      }
      
      this->eeComTask_->A().setFromTriplets(eeComTripletList_A.begin(), eeComTripletList_A.end());
      this->eeComTask_->C().setFromTriplets(eeComTripletList_C.begin(), eeComTripletList_C.end());
      this->eeComTask_->wa() = Eigen::VectorXd::Ones(6 * gaitParam.eeName.size() + 3)* 1e-3;
      this->eeComTask_->dl() = Eigen::VectorXd::Zero(6 + actRobotTqc->numJoints());
      this->eeComTask_->du() = Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * 3.0;
      this->eeComTask_->wc() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints()) * 1e-6;
      this->eeComTask_->w() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints()) * 1e-6;
      this->eeComTask_->toSolve() = true;
      this->eeComTask_->settings().check_termination = 15; // default 25. 高速化
      this->eeComTask_->settings().verbose = 0;

      // rootTask
      // Jacobianもgenerate frameなのでdJ * dq = 0
      {
	// rootTaskを作る
	this->rootTask_->A() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6,6 + actRobotTqc->numJoints());
	std::vector<Eigen::Triplet<double> > tripletList_A;
	for (int i=0;i<6;i++){
	  tripletList_A.push_back(Eigen::Triplet<double>(i,i,1.0));
	}
	this->rootTask_->b() = cnoid::Vector6::Zero();
	for (int i=0;i<3;i++){
	  this->rootTask_->b()[i] = this->root_K[i] * (gaitParam.genRobot->rootLink()->p()[i] - gaitParam.actRobot->rootLink()->p()[i]);
	}
	for (int i=0;i<3;i++){
	  this->rootTask_->b()[i+3] = this->root_K[i+3] * cnoid::rpyFromRot(gaitParam.genRobot->rootLink()->R() * gaitParam.actRobot->rootLink()->R().transpose())[i];
	}
	
	this->rootTask_->wa() = cnoid::Vector6::Ones() * 1e-3;
	this->rootTask_->C() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 + actRobotTqc->numJoints(),6 + actRobotTqc->numJoints());
	std::vector<Eigen::Triplet<double> > tripletList_C;
	for (int i=0;i<6;i++){
	  tripletList_C.push_back(Eigen::Triplet<double>(i,i,1.0));
	}
	this->rootTask_->dl() = - Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * 1e+1;
	this->rootTask_->du() = Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * 1e+1;
	this->rootTask_->wc() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints()) * 1e-3;
	this->rootTask_->w() = cnoid::VectorXd::Ones(6 + actRobotTqc->numJoints()) * 1e-6;
	this->rootTask_->toSolve() = true;
	this->rootTask_->settings().check_termination = 15; // default 25. 高速化
	this->rootTask_->settings().verbose = 0;
      }      

      // ddqを計算
      std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
      tasks.push_back(this->eeComTask_);
      tasks.push_back(this->rootTask_);
	
      cnoid::VectorX result;
      if(!prioritized_qp_base::solve(tasks,
				     result,
				     0 // debuglevel
				     )){
	std::cerr << "fail" << std::endl; // TODO
      }else {
	// ddqを代入
	actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	actRobotTqc->rootLink()->dv() = result.block(0,0,3,1);//  << result[0], result[1], result[2];
	actRobotTqc->rootLink()->dw() = result.block(3,0,3,1);// << result[3], result[4], result[5];
	for(int i=0;i<actRobotTqc->numJoints();i++){
	  actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	  // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	  // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	  actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	  actRobotTqc->joint(i)->ddq() = result[6+i];
	}
	cnoid::calcInverseDynamics(actRobotTqc->rootLink()); // actRobotTqc->joint()->u()に書き込まれる
	for(int i=0;i<actRobotTqc->numJoints();i++){
	  tau_ee[i] = actRobotTqc->joint(i)->u();
	  actRobotTqc->joint(i)->u() = 0.0;
	}
      }      
    }

    // 最終的な出力トルクを代入
    for(int i=0;i<actRobotTqc->numJoints();i++){
      actRobotTqc->joint(i)->u() = tau_g[i] + tau_lip[i] + tau_ee[i];
    }

    /*
    std::cerr << "tau_g" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_g[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau_lip" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_lip[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau_ee" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_ee[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << actRobotTqc->joint(i)->u() << " ";
    }
    std::cerr << std::endl;*/
    
    // Gain
    {
      for(int i=0;i<NUM_LEGS;i++){
	cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	if(gaitParam.isManualControlMode[i].getGoal() == 0.0) { // Manual Control off
	  if(gaitParam.footstepNodesList[0].isSupportPhase[i]){
	    double transitionTime = std::max(this->landing2SupportTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->supportPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->supportPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->supportDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->supportDgain[i][j], transitionTime);
	    }
	  }else if(gaitParam.swingState[i] == GaitParam::DOWN_PHASE) {
	    double transitionTime = std::max(this->swing2LandingTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->landingPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->landingPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->landingDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->landingDgain[i][j], transitionTime);
	    }
	  }else{
	    double transitionTime = std::max(this->support2SwingTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->swingPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->swingPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->swingDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->swingDgain[i][j], transitionTime);
	    }
	  }
	}else{ // Manual Control on
	  double transitionTime = std::max(gaitParam.isManualControlMode[i].remain_time(), dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	  for(int j=0;j<jointPath.numJoints();j++){
	    if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != 100.0) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(100.0, transitionTime);
	    if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != 100.0) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(100.0, transitionTime);
	  }
	}
      }

      // arm
      for(int i=NUM_LEGS;i<NUM_LEGS+2;i++){
	cnoid::JointPath jointPath(actRobotTqc->link("CHEST_JOINT2"), actRobotTqc->link(gaitParam.eeParentLink[i]));
	double arm_gain = 0.0;
	for(int j=0;j<jointPath.numJoints();j++){
	  if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != arm_gain) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(arm_gain, 3.0);
	  if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != arm_gain) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(arm_gain, 3.0);
	}
      }
    }
  
    for(int i=0;i<gaitParam.genRobot->numJoints();i++){
      o_stServoPGainPercentage[i].interpolate(dt);
      o_stServoDGainPercentage[i].interpolate(dt);
    }
  } // useActState

  for(int i=0;i<gaitParam.genRobot->numJoints();i++){
    prev_q[i] = gaitParam.genRobot->joint(i)->q();
    prev_dq[i] = gaitParam.genRobot->joint(i)->dq();
  }

  for(int i=0;i<gaitParam.eeName.size();i++){
    eeTargetPosedd[i] = eeTargetPosed[i];
    eeTargetPosed[i] = gaitParam.abcEETargetPose[i];
  }
  
  return true;
}

