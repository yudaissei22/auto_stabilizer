#include "LegCoordsGenerator.h"
#include "FootGuidedController.h"
#include "MathUtil.h"

namespace legcoordsgenerator{

  void calcLegCoords(GaitParam& gaitParam, double dt){
    // swing期は、remainTime - supportTime - delayTimeOffset後にdstCoordsに到達するようなantececdent軌道を生成し(genCoords.getGoal()の値)、その軌道にdelayTimeOffset遅れで滑らかに追従するような軌道(genCoords.value()の値)を生成する.
    //   rectangle以外の軌道タイプや跳躍についてはひとまず考えない TODO
    //   srcCoordsとdstCoordsのsrcCoordsとdstCoordsの高い方よりもさらにstepHeightだけ高い位置(heightとおく)に上げるようなrectangle軌道を生成する
    //     dstに到達するまでの時間が(height-dstCoords)/touchVel以下の場合、dstまで直線で移動する
    //     そうで無い場合、今の位置がheight - eps(ひとまず0)よりも低ければ、heightの高さまで上げてからdstCoordsの上空(XYと回転はdstCoordsと同じ)まで直線で移動する軌道を,残りの自由に使える時間で実行する
    //                     今の位置がheight - eps(ひとまず0)以上であれば、dstCoordsの上空(XYと回転はdstCoordsと同じ)まで直線で移動する軌道を,残りの自由に使える時間で実行する
    //                     両者はチャタリングしそうだが、delayTimeOffset遅れで滑らかに追従するので大丈夫
    // support期は、現FootStepNodesの終了時にdstCoordsに到達するような軌道を線形補間によって生成する.

    for(int i=0;i<2;i++){
      if(gaitParam.footstepNodesList[0].remainTime <= gaitParam.footstepNodesList[0].supportTime[i]) { // 支持脚
        cnoid::Position nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{gaitParam.genCoords[i].value(),gaitParam.footstepNodesList[0].dstCoords[i]},
                                             std::vector<double>{std::max(0.0,gaitParam.footstepNodesList[0].remainTime - dt), dt}); // このfootstepNode終了時にdstCoordsに行くように線形補間
        gaitParam.genCoords[i].reset(nextCoords);
      }else{ // 遊脚
        double swingTime = gaitParam.footstepNodesList[0].remainTime - gaitParam.footstepNodesList[0].supportTime[i];
        double height = std::max(gaitParam.srcCoords[i].translation()[2],gaitParam.footstepNodesList[0].dstCoords[i].translation()[2]) + gaitParam.footstepNodesList[0].stepHeight[i]; // 足上げ高さ. generate frame
        cnoid::Position srcCoords = gaitParam.srcCoords[i];
        cnoid::Position dstCoords = gaitParam.footstepNodesList[0].dstCoords[i];
        cnoid::Position antecedentCoords = gaitParam.genCoords[i].getGoal(); // 今のantecedent軌道の位置
        if(swingTime <= gaitParam.delayTimeOffset) {
          gaitParam.genCoords[i].setGoal(dstCoords, swingTime);
          gaitParam.genCoords[i].interpolate(dt);
        }else{
          double antecedentRemainTime = swingTime - gaitParam.delayTimeOffset; // antecedent軌道はこの時間後にdstCoordsにつく
          double touchDownTime = (antecedentCoords.translation()[2] - dstCoords.translation()[2]) / gaitParam.touchVel; // 地面につくのに要する時間
          cnoid::Position nextCoords;
          if(antecedentRemainTime <= touchDownTime){ // 地面に直線的につきにいく
             nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{antecedentCoords,dstCoords},
                                                  std::vector<double>{std::max(0.0,antecedentRemainTime - dt), dt}); // このfootstepNode終了時にdstCoordsに行くように線形補間
          }else{
            double preTouchDownRemainTime = antecedentRemainTime - touchDownTime; // dstCoords の上空につくまでの時間
            cnoid::Position preTouchCoords = dstCoords; // dstCoords の上空
            preTouchCoords.translation()[2] = height;
            if(antecedentCoords.translation()[2] >= height){ // dstCoords の上空まで直線で移動する
              nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{antecedentCoords,preTouchCoords},
                                                   std::vector<double>{std::max(0.0,preTouchDownRemainTime - dt), dt}); // preTouchDownRemainTime後にpreTouchCoordsに行くように線形補間
            }else{ // heightの高さまで上げてからdstCoords の上空まで直線で移動する. 回転は2つの区間の間一定の速度で動き続ける
              cnoid::Vector3 viaPos = antecedentCoords.translation(); // heightの高さまで上げた位置
              viaPos[2] = height;
              double viaTime = preTouchDownRemainTime * (viaPos - antecedentCoords.translation()).norm() / ((viaPos - antecedentCoords.translation()).norm() + (preTouchCoords.translation() - viaPos).norm()); // viaPosにつくまでの時間
              nextCoords.translation() = mathutil::calcMidPos(std::vector<cnoid::Vector3>{antecedentCoords.translation(),viaPos},
                                                              std::vector<double>{std::max(0.0,viaTime - dt), dt}); // viaTime後にviaPosに行くように線形補間
              nextCoords.linear() =  mathutil::calcMidRot(std::vector<cnoid::Matrix3>{antecedentCoords.linear(),preTouchCoords.linear()},
                                                          std::vector<double>{std::max(0.0,preTouchDownRemainTime - dt), dt}); // preTouchDownRemainTime後にpreTouchCoordsの傾きになるように線形補間
            }
          }
          gaitParam.genCoords[i].setGoal(nextCoords, gaitParam.delayTimeOffset);
          gaitParam.genCoords[i].interpolate(dt);
        }
      }
    }

    // refZmp TODO ここから

    gaitParam.footstepNodesList[0].remainTime -= dt;
    if(gaitParam.footstepNodesList[0].remainTime <= 0.0 && gaitParam.footstepNodesList[0].remainTime.size() > 1){
      for(int i=0;i<2;i++){
        gaitParam.srcCoords[i] = gaitParam.genCoords[i].value();
      }
      gaitParam.footstepNodesList.erase(gaitParam.footstepNodesList.begin()); // vectorではなくlistにするべき?
    }
  }

  void calcCOMCoords(const GaitParam& gaitParam, double dt, double g, double mass, cnoid::Vector3& genNextCog, cnoid::Vector3& genNextCogVel){
    std::vector<footguidedcontroller::LinearTrajectory<cnoid::Vector3> > refZmpTraj;
    {
      std::vector<cnoid::Position> stepCoords{gaitParam.genCoords[0].value(),gaitParam.genCoords[1].value()}; // for文中の現在の脚の位置
      cnoid::Vector3 refZmp = gaitParam.refZmp; // for文中の現在のrefzmp
      for(int i=0;i<gaitParam.footstepNodesList.size();i++){
        double remainTime = gaitParam.footstepNodesList[i].remainTime;

        // 跳躍についてはひとまず考えない TODO
        // if(gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[0] &&
        //    gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[1])

        if(remainTime > gaitParam.footstepNodesList[i].supportTime[0] && remainTime <= gaitParam.footstepNodesList[i].supportTime[1]){// 右脚がswing. refzmpは左脚の位置
          double swingTime = remainTime - gaitParam.footstepNodesList[i].supportTime[0];
          cnoid::Position llegStartCoords = stepCoords[1];
          cnoid::Position llegGoalCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{stepCoords[1],gaitParam.footstepNodesList[i].dstCoords[1]},std::vector<double>{swingTime, remainTime-swingTime}); // このfootstepNode終了時にdstCoordsに行くように線形補間
          cnoid::Vector3 zmpStartPos = llegStartCoords.translation() + llegStartCoords.linear()*gaitParam.copOffset[1];
          cnoid::Vector3 zmpGoalPos = llegGoalCoords.translation() + llegGoalCoords.linear()*gaitParam.copOffset[1];
          refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,swingTime));
          stepCoords[0] = gaitParam.footstepNodesList[i].dstCoords[0];
          stepCoords[1] = llegGoalCoords;
          refZmp = zmpGoalPos;
          remainTime -= swingTime;
        }else if(gaitParam.footstepNodesList[i].remainTime <= gaitParam.footstepNodesList[i].supportTime[0] && gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[1]){ // 左脚がswing. refzmpは右脚の位置
          double swingTime = remainTime - gaitParam.footstepNodesList[i].supportTime[1];
          cnoid::Position rlegStartCoords = stepCoords[0];
          cnoid::Position rlegGoalCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{stepCoords[0],gaitParam.footstepNodesList[i].dstCoords[0]},std::vector<double>{swingTime, remainTime-swingTime}); // このfootstepNode終了時にdstCoordsに行くように線形補間
          cnoid::Vector3 zmpStartPos = rlegStartCoords.translation() + rlegStartCoords.linear()*gaitParam.copOffset[0];
          cnoid::Vector3 zmpGoalPos = rlegGoalCoords.translation() + rlegGoalCoords.linear()*gaitParam.copOffset[0];
          refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,swingTime));
          stepCoords[0] = rlegGoalCoords;
          stepCoords[1] = gaitParam.footstepNodesList[i].dstCoords[1];
          refZmp = zmpGoalPos;
          remainTime -= swingTime;
        }
        if(remainTime <= gaitParam.footstepNodesList[i].supportTime[0] && remainTime <= gaitParam.footstepNodesList[i].supportTime[1]){ // 両脚がsupport. refzmpは一つ前の区間と一つ後の区間を線形につなぐ
          cnoid::Position rlegGoalCoords = gaitParam.footstepNodesList[i].dstCoords[0];
          cnoid::Position llegGoalCoords = gaitParam.footstepNodesList[i].dstCoords[1];
          cnoid::Vector3 rlegCOP = rlegGoalCoords.translation() + rlegGoalCoords.linear()*gaitParam.copOffset[0];
          cnoid::Vector3 llegCOP = llegGoalCoords.translation() + llegGoalCoords.linear()*gaitParam.copOffset[1];
          cnoid::Vector3 zmpStartPos = refZmp;
          cnoid::Vector3 zmpGoalPos;
          if(i==gaitParam.footstepNodesList.size()-1 || //末尾. 以降は末尾の状態がずっと続くとして扱うので、refzmpは両足の中心
             gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[0] && gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[1] // 次も両脚支持. refzmpは両足の中心
             ){
            zmpGoalPos = 0.5 * rlegCOP + 0.5 * llegCOP;
          }else if(gaitParam.footstepNodesList[i+1].remainTime > gaitParam.footstepNodesList[i+1].supportTime[0] && gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[1]){ // 次は右脚がswing
            zmpGoalPos = llegCOP;
          }else if(gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[0] && gaitParam.footstepNodesList[i+1].remainTime > gaitParam.footstepNodesList[i+1].supportTime[1]){ // 次は左脚がswing
            zmpGoalPos = rlegCOP;
          }else{// 跳躍についてはひとまず考えない TODO
          }
          refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,gaitParam.footstepNodesList[i].remainTime));
          stepCoords[0] = rlegGoalCoords;
          stepCoords[1] = llegGoalCoords;
          refZmp = zmpGoalPos;
          remainTime = 0.0;
        }
      }
    }

    double w = std::sqrt(g/gaitParam.dz);
    cnoid::Vector3 l = cnoid::Vector3::Zero();
    l[2] = gaitParam.dz;
    cnoid::Vector3 genDCM = gaitParam.genCog + w * gaitParam.genCogVel;
    cnoid::Vector3 genZmp = footguidedcontroller::calcFootGuidedControl(w,l,genDCM,refZmpTraj);
    // check zmp in polygon TODO
    cnoid::Vector3 genNextForce;
    footguidedcontroller::updateState(w,l,gaitParam.genCog,gaitParam.genCogVel,genZmp,mass,dt,
                                      genNextCog, genNextCogVel, genNextForce);
    return;
  }

  bool calcNextCoords(GaitParam& gaitParam, double dt, double g, double mass){
    calcLegCoords(gaitParam,dt);
    cnoid::Vector3 genNextCog,genNextCogVel;
    calcCOMCoords(gaitParam, dt, g, mass,genNextCog,genNextCogVel);
    gaitParam.genCog = genNextCog;
    gaitParam.genCogVel = genNextCogVel;
    return true;
  }
}
